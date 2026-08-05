#!/usr/bin/env python3
import argparse
import asyncio
import json
import os
import re
import pandas as pd
from playwright.async_api import async_playwright, TimeoutError as PWTimeoutError

URL = "https://clingen.igib.res.in/indigen/index.php"
INPUT_SELECTOR = 'input[placeholder="Variant / Gene / dbSNP ID"]'

CHR_ORDER = [
    "Chr", "Position", "Ref Allele", "Alt Allele", "Gene",
    "Exonic Function", "Allele Count", "Allele Frequency",
    "Allele Number", "Homozygous", "Heterozygous",
]  # preferred first columns if these exact names show up; anything else
   # the API returns is kept too, appended after, so nothing is lost

# One generous wait per gene, instead of several short retries.
# Big genes (lots of variants) just need more time server-side -
# retrying restarts the query rather than giving it more time.
MAX_TOTAL_WAIT_S = 15 * 60   # 15 minute ceiling per gene, as requested
FIRST_TIMEOUT_S = 45         # quick first check - most genes finish fast
FALLBACK_TIMEOUT_S = MAX_TOTAL_WAIT_S - FIRST_TIMEOUT_S  # ~14min patient wait

# route.fetch() itself must be bounded too, or a genuinely stuck server-side
# query leaves an orphaned fetch running in the background forever (even
# after we've given up on this gene client-side), which backs up
# Playwright's driver process and corrupts its internal timeout tracking -
# that's what caused the "TimeoutNegativeWarning" / negative timeout errors.
ROUTE_FETCH_TIMEOUT_MS = (MAX_TOTAL_WAIT_S + 20) * 1000

ROW_COUNT_REFRESH_THRESHOLD = 1000  # refresh page after big result sets
PERIODIC_RESTART_EVERY = 25         # full browser relaunch every N genes,
                                     # proactively, to avoid slow memory
                                     # creep before it causes a crash

CHROME_LAUNCH_ARGS = [
    # headless Chromium on Linux frequently crashes rendering large DOMs
    # because /dev/shm defaults to a small size - this is the standard fix
    "--disable-dev-shm-usage",
    "--disable-gpu",
]


def safe_filename(gene):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", gene)


def parse_rows(obj):
    """
    Return each variant record AS-IS from the API (minus the internal
    Mongo _id), instead of mapping to a fixed, guessed set of column
    names. Different key names (or missing keys for different variant
    types - SNV vs indel, coding vs non-coding) are preserved exactly
    as the site sends them, so no field silently goes missing.
    """
    rows = []
    for x in obj.get("mydata", []):
        row = dict(x)
        row.pop("_id", None)
        rows.append(row)
    return rows


class GeneSearcher:
    """
    Registers ONE persistent route handler for data.php on the page and
    reuses it across all genes. Each search creates a fresh Future; the
    route handler only resolves that Future if the response actually
    corresponds to the gene we just asked for (checked via URL/post data),
    so a late/stale response from a previous slow gene can never get
    attributed to the wrong gene.

    Body is read via route.fetch() -> APIResponse.text(), which goes
    through Playwright's own network layer, NOT Chrome's DevTools
    inspector cache - this avoids the "evicted from inspector cache"
    error entirely, since that cache is never involved.
    """

    def __init__(self, page):
        self.page = page
        self.current_gene = None
        self.future = None
        self.generation = 0  # bumped every new search; lets us detect and
                              # abort requests we've since abandoned instead
                              # of letting them deliver a stale/huge payload
                              # to the page at an unpredictable later time

    async def start(self):
        await self.page.route("**/data.php*", self._handle_route)

    async def stop(self):
        try:
            await self.page.unroute("**/data.php*", self._handle_route)
        except Exception:
            pass

    async def rebind(self, new_page):
        """Point this searcher at a fresh page (new tab / clean DOM)."""
        self.page = new_page
        self.current_gene = None
        self.future = None
        self.generation += 1
        await self.start()

    async def _handle_route(self, route, request):
        my_generation = self.generation

        if my_generation != self.generation:
            try:
                await route.abort()
            except Exception:
                pass
            return

        response = None
        try:
            response = await route.fetch(timeout=ROUTE_FETCH_TIMEOUT_MS)
            body = await response.text()
        except Exception as e:
            body = None
            fetch_error = e
        else:
            fetch_error = None

        # If we've since moved on (new gene, retry, or abandoned this one)
        # while that fetch was running, don't deliver it to the page now -
        # a late, possibly huge response landing unexpectedly mid-interaction
        # is what was freezing the page and corrupting timeouts. Abort it
        # cleanly instead.
        if my_generation != self.generation:
            try:
                await route.abort()
            except Exception:
                pass
            return

        try:
            if response is not None:
                await route.fulfill(response=response)
            else:
                await route.continue_()
        except Exception:
            pass

        if self.future is None or self.future.done():
            return

        gene = (self.current_gene or "").upper()
        haystack = (request.url or "")
        try:
            post_data = request.post_data or ""
        except Exception:
            post_data = ""
        haystack += post_data
        haystack = haystack.upper()

        if gene and gene not in haystack:
            return  # response for a different (likely stale) query - ignore

        if fetch_error is not None:
            if not self.future.done():
                self.future.set_exception(fetch_error)
            return

        if not self.future.done():
            self.future.set_result(body)

    async def search(self, gene):
        self.current_gene = gene
        self.future = asyncio.get_event_loop().create_future()
        self.generation += 1

        try:
            locator = self.page.locator(INPUT_SELECTOR)
            await locator.scroll_into_view_if_needed(timeout=30000)
            await locator.click(timeout=30000)
            await locator.fill("", timeout=15000)
            await locator.fill(gene, timeout=15000)
            await self.page.keyboard.press("Enter")
        except Exception as e:
            self.current_gene = None
            raise  # let caller decide if this means the browser died

        # Submit ONCE, then just wait longer on the same in-flight request -
        # never resubmit the same gene. Use asyncio.wait() (not wait_for),
        # because wait_for CANCELS the awaited future the moment it times
        # out - reusing that same (now-cancelled) future for a second,
        # longer wait would raise CancelledError immediately instead of
        # actually waiting. asyncio.wait() leaves an unfinished future
        # alone on timeout, so we can genuinely keep waiting on it.
        done, _ = await asyncio.wait({self.future}, timeout=FIRST_TIMEOUT_S)
        if self.future not in done:
            print(f"{gene:20s}  still running, giving it up to {MAX_TOTAL_WAIT_S}s total...")
            done, _ = await asyncio.wait(
                {self.future}, timeout=MAX_TOTAL_WAIT_S - FIRST_TIMEOUT_S
            )
            if self.future not in done:
                self.current_gene = None
                return None

        try:
            body = self.future.result()
        except Exception as e:
            print(f"    fetch error for {gene}: {e}")
            self.current_gene = None
            return None

        self.current_gene = None
        if not body or not body.strip():
            return None
        try:
            obj = json.loads(body)
        except Exception as e:
            print(f"    bad JSON for {gene}: {e}")
            return None
        return parse_rows(obj)


async def search_gene(searcher, gene):
    rows = await searcher.search(gene)
    if rows is not None:
        print(f"{gene:20s} {len(rows):6d}")
    else:
        print(f"{gene:20s}  FAILED (no response within {MAX_TOTAL_WAIT_S}s)")
    return rows


def is_browser_dead_error(e):
    msg = str(e).lower()
    return "has been closed" in msg or "target closed" in msg or "browser has disconnected" in msg


async def launch_browser(p):
    browser = await p.chromium.launch(headless=True, args=CHROME_LAUNCH_ARGS)
    page = await browser.new_page()
    await page.goto(URL, wait_until="networkidle")
    return browser, page


def save(all_rows, outpath):
    if not all_rows:
        pd.DataFrame().to_csv(outpath, sep="\t", index=False)
        return pd.DataFrame()

    df = pd.DataFrame(all_rows)  # columns = union of every key seen, in
                                  # order of first appearance - nothing
                                  # the API returned gets dropped

    preferred = [c for c in CHR_ORDER if c in df.columns]
    rest = [c for c in df.columns if c not in preferred]
    df = df[preferred + rest]

    try:
        df.drop_duplicates(inplace=True)
    except TypeError:
        dedup_key = df.astype(str).apply(tuple, axis=1)
        df = df[~dedup_key.duplicated()]

    df.to_csv(outpath, sep="\t", index=False)
    return df


def rebuild_output_from_cache(genes, cache_dir, outpath):
    """
    Rebuild the full output file from whatever is cached on disk,
    covering this run AND any previous runs. This is what makes the
    script resumable: re-running the exact same command will only
    fetch genes that are still missing, and the final file always
    reflects everything ever successfully cached.
    """
    all_rows = []
    missing = []
    for gene in genes:
        p = os.path.join(cache_dir, safe_filename(gene) + ".json")
        if os.path.exists(p):
            try:
                with open(p) as f:
                    all_rows.extend(json.load(f))
            except Exception as e:
                print(f"    couldn't read cache for {gene}, will re-fetch: {e}")
                missing.append(gene)
        else:
            missing.append(gene)
    df = save(all_rows, outpath)
    return df, missing


async def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-gl", required=True)
    parser.add_argument("-o", required=True)
    args = parser.parse_args()

    with open(args.gl) as f:
        genes = [x.strip() for x in f if x.strip()]

    cache_dir = args.o + ".cache"
    os.makedirs(cache_dir, exist_ok=True)

    already_cached = [g for g in genes if os.path.exists(os.path.join(cache_dir, safe_filename(g) + ".json"))]
    todo = [g for g in genes if g not in already_cached]
    if already_cached:
        print(f"{len(already_cached)} gene(s) already cached from a previous run, skipping them: "
              f"{', '.join(already_cached[:10])}{' ...' if len(already_cached) > 10 else ''}")
    print(f"{len(todo)} gene(s) to fetch this run.")
    print()

    async with async_playwright() as p:
        browser, page = await launch_browser(p)
        searcher = GeneSearcher(page)
        await searcher.start()

        for i, gene in enumerate(todo, 1):
            rows = None
            for browser_attempt in range(2):  # normal try, then one full relaunch+retry
                try:
                    rows = await search_gene(searcher, gene)
                    break
                except Exception as e:
                    if is_browser_dead_error(e) and browser_attempt == 0:
                        print(f"    browser died on {gene}, relaunching...")
                        try:
                            await browser.close()
                        except Exception:
                            pass
                        try:
                            browser, page = await launch_browser(p)
                            await searcher.rebind(page)
                            continue  # retry this gene once on the fresh browser
                        except Exception as e2:
                            print(f"    relaunch failed: {e2}")
                            rows = None
                            break
                    else:
                        print(f"{gene:20s}  unexpected error: {e}")
                        rows = None
                        break

            if rows is not None:
                cache_path = os.path.join(cache_dir, safe_filename(gene) + ".json")
                with open(cache_path, "w") as f:
                    json.dump(rows, f)

            needs_refresh = (rows is None) or (len(rows) > ROW_COUNT_REFRESH_THRESHOLD) \
                or (i % PERIODIC_RESTART_EVERY == 0)

            if needs_refresh:
                print(f"    refreshing browser (clean state) before next gene...")
                try:
                    await browser.close()
                except Exception:
                    pass
                try:
                    browser, page = await launch_browser(p)
                    await searcher.rebind(page)
                except Exception as e:
                    print(f"    couldn't refresh browser: {e}")
            else:
                await page.wait_for_timeout(1500)

            if i % 10 == 0:
                rebuild_output_from_cache(genes, cache_dir, args.o)

        await searcher.stop()
        try:
            await browser.close()
        except Exception:
            pass

    df, missing = rebuild_output_from_cache(genes, cache_dir, args.o)

    print()
    print("Variants:", len(df))
    if missing:
        print(f"Still missing ({len(missing)}) - just re-run the same command to retry only these:")
        print(", ".join(missing))
    else:
        print("All genes fetched successfully.")


if __name__ == "__main__":
    asyncio.run(main())
