#!/usr/bin/env python3

import argparse
import csv
import os
import time

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options

from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.chrome.service import Service

from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

from selenium.common.exceptions import TimeoutException, NoSuchElementException


BASE_URL = "https://clingen.igib.res.in/indigen/index.php"


# -------------------------------------------------
# Arguments
# -------------------------------------------------
def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument("-gl", "--genefile", required=True)
    parser.add_argument("-o", "--output", required=True)

    return parser.parse_args()


# -------------------------------------------------
# Read genes
# -------------------------------------------------
def read_gene_list(file):

    genes = []

    with open(file) as f:
        for line in f:
            g = line.strip()
            if g:
                genes.append(g)

    return genes


# -------------------------------------------------
# Setup browser
# -------------------------------------------------
def setup_driver():

    chrome_binary = "/usr/bin/google-chrome"

    if not os.path.exists(chrome_binary):
        raise RuntimeError("Google Chrome not found")

    print("Using browser:", chrome_binary)

    options = Options()
    options.binary_location = chrome_binary

    options.add_argument("--headless=new")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--window-size=1920,1080")

    service = Service(ChromeDriverManager().install())

    driver = webdriver.Chrome(
    service=service,
    options=options,
    keep_alive=True,
    )
    driver.set_page_load_timeout(600)
    driver.set_script_timeout(600)

    return driver

# -------------------------------------------------
# Extract variant data
# -------------------------------------------------
def extract_variant_data(driver, gene):

    body = driver.find_element(By.TAG_NAME, "body").text

    lines = body.split("\n")

    data = {}

    i = 0

    while i < len(lines):

        line = lines[i].strip()

        if ":" in line:

            key, val = line.split(":", 1)

            key = key.strip()
            val = val.strip()

            if val == "" and i + 1 < len(lines):
                val = lines[i + 1].strip()
                i += 1

            data[key] = val

        i += 1

    return {
        "Chr": data.get("Chr", ""),
        "Position": data.get("Position", ""),
        "Ref Allele": data.get("Ref Allele", ""),
        "Alt Allele": data.get("Alt Allele", ""),
        "Gene": gene,
        "Exonic Function": data.get("Exonic Function", ""),
        "Allele Count": data.get("Allele Count", ""),
        "Allele Frequency": data.get("Allele Frequency", ""),
        "Allele Number": data.get("Allele Number", ""),
        "Homozygous": data.get("Homozygous", ""),
        "Heterozygous": data.get("Heterozygous", "")
    }


# -------------------------------------------------
# Collect variant links
# -------------------------------------------------
def collect_variant_links(driver):

    wait = WebDriverWait(driver, 300)

    links = set()
    last_count = 0

    while True:

        wait.until(
            EC.presence_of_element_located(
                (By.XPATH, "//a[contains(@href,'showdata.php')]")
            )
        )

        elements = driver.find_elements(
            By.XPATH, "//a[contains(@href,'showdata.php')]"
        )

        for e in elements:
            links.add(e.get_attribute("href"))

        print("Collected so far:", len(links))

        if len(links) == last_count:
            print("No new variants found. Stopping pagination.")
            break

        last_count = len(links)

        try:

            next_button = driver.find_element(
                By.XPATH, "//a[contains(text(),'Next')]"
            )

            if "disabled" in next_button.get_attribute("class"):
                break

            driver.execute_script("arguments[0].click();", next_button)

            time.sleep(2)

        except NoSuchElementException:
            break

    return list(links)


# -------------------------------------------------
# Process gene
# -------------------------------------------------
def process_gene(driver, gene, writer=None):

    results = []

    driver.get(BASE_URL)

    wait = WebDriverWait(driver, 300)

    try:

        search_box = wait.until(
            EC.presence_of_element_located((By.XPATH, "//input[@type='text']"))
        )

        search_box.clear()
        search_box.send_keys(gene)
        search_box.send_keys(Keys.ENTER)

        wait.until(
            EC.presence_of_element_located(
                (By.XPATH, "//a[contains(@href,'showdata.php')]")
            )
        )

        try:
            urls = collect_variant_links(driver)
        except Exception as e:
            print(f"Failed for {gene}: {e}")
            return []

        print("Total variants found:", len(urls))

        for i, url in enumerate(urls):

            print(f"Downloading variant {i+1}/{len(urls)}")

            driver.get(url)

            wait.until(
                EC.presence_of_element_located(
                    (By.XPATH, "//*[contains(text(),'IndiGen')]")
                )
            )

            data = extract_variant_data(driver, gene)
            results.append(data)

            if writer and (i + 1) % 50 == 0:
                writer.writerows(results)
                results = []

    except TimeoutException:
        print("No variants found for", gene)

    return results


# -------------------------------------------------
# Main
# -------------------------------------------------
def main():

    args = parse_arguments()

    genes = read_gene_list(args.genefile)

    driver = setup_driver()

    fields = [
        "Chr",
        "Position",
        "Ref Allele",
        "Alt Allele",
        "Gene",
        "Exonic Function",
        "Allele Count",
        "Allele Frequency",
        "Allele Number",
        "Homozygous",
        "Heterozygous"
    ]

    total = 0

    with open(args.output, "w", newline="") as f:

        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()

        try:

            for gene in genes:

                print("\nProcessing", gene)

                rows = process_gene(driver, gene, writer)

                writer.writerows(rows)

                total += len(rows)

        finally:
            driver.quit()

    print("\nFinished")
    print("Total variants:", total)
    print("Output written:", args.output)


if __name__ == "__main__":
    main()
