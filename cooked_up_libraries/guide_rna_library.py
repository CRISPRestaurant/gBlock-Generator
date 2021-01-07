import os
import time
import shutil
from datetime import date, datetime
import csv
import pandas as pd
import sys

from selenium import webdriver
from selenium.webdriver.support.ui import Select, WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException, NoSuchElementException
from selenium.webdriver.common.by import By

# Uses the CHOPCHOP API to obtain guideRNA candidates for an
# organism and a gene from the offline runnable CHOPCHOP Python API

# WARNING: This is complete hell to setup so use this as a resource (https://github.com/CRISPRestaurant/gBlock-Generator/blob/master/Instructions/guideRNA%20Generation%20Python%20Tutorial.md)
# if you want to be able to run it. This setup only works on Linux.
def obtainGuideRNAsOfflineCHOPCHOP(organism_name, target_gene, guide_rna_results_storage_folder, chopchop_dir):
    # Removes all previous off target data downloaded into the off_target_site_destination_folder
    organism_folder = "%s/%s" % (guide_rna_results_storage_folder, organism_name)
    gene_folder = "%s/%s/%s" % (guide_rna_results_storage_folder, organism_name, target_gene)
    off_target_folder = "%s/%s/%s/Off_Targets" % (guide_rna_results_storage_folder, organism_name, target_gene)
    if os.path.isdir(organism_folder):
        if os.path.isdir(gene_folder):
            shutil.rmtree(gene_folder)

        os.makedirs(off_target_folder)
    else:
        os.makedirs(off_target_folder)

    # END

    command = "%s/chopchop.py -G %s -o %s -Target %s" % (chopchop_dir, organism_name, off_target_folder, target_gene)
    command_stream = os.popen(command)
    command_output = command_stream.read()

    command_output_split = command_output.split("\n")

    for i in range(len(command_output_split)):
        command_output_split[i] = command_output_split[i].split("\t")

    command_output_split = command_output_split[0: len(command_output_split) - 1]

    # Generates output summary of run
    guide_rna_output_summary_file = open("%s/Knockout Output Summary.txt" % (gene_folder), "w")
    guide_rna_output_summary_file.write("Date Run: %s" % (date.today()))
    guide_rna_output_summary_file.close()

    # Saves output, a 2D array, into csv format for future retrieval
    column_headers = command_output_split[0]
    non_column_headers = command_output_split[1:]

    with open("%s/guideRNAs.csv" % (gene_folder), "w") as f:
        write = csv.writer(f)

        write.writerow(column_headers)
        write.writerows(non_column_headers)

    return command_output_split

# Obtains guideRNAs using Selenium and CHOPCHOP website
def obtainGuideRNAsOnlineCHOPCHOP(organism_name, target_gene, guide_rna_results_storage_folder, time_allotment = 120, selenium_driver = "chromedriver/chromedriver"):
    # Removes all previous off target data downloaded into the off_target_site_destination_folder
    organism_folder = "%s/%s" % (guide_rna_results_storage_folder, organism_name)
    gene_folder = "%s/%s/%s" % (guide_rna_results_storage_folder, organism_name, target_gene)
    off_target_folder = "%s/%s/%s/Off_Targets" % (guide_rna_results_storage_folder, organism_name, target_gene)
    if os.path.isdir(organism_folder):
        if os.path.isdir(gene_folder):
            shutil.rmtree(gene_folder)

        os.makedirs(off_target_folder)
    else:
        os.makedirs(off_target_folder)

    # END

    output = [['Rank', 'Target sequence', 'Genomic location', 'Strand', 'GC content (%)', 'Self-complementarity', 'MM0', 'MM1', 'MM2', 'MM3', 'Efficiency']]

    driver = webdriver.Chrome(executable_path = selenium_driver)

    driver.get("https://chopchop.cbu.uib.no/")

    # Inputs organism name by finding it among selection and inputs gene name into respective input field
    target_gene_element = driver.find_element_by_xpath("//*[@id=\"geneInput\"]")
    target_species_element = Select(driver.find_element_by_xpath("//*[@id=\"genomeSelect\"]"))

    target_species_options = target_species_element.options

    target_species_value_matched_with_input = ""
    for target_species_option in target_species_options:
        target_species_candidate = target_species_option.text
        simplified_names = target_species_candidate[target_species_candidate.find("(") + 1: target_species_candidate.find(")")]

        simplified_names_list = [simplified_names]

        if ", " in simplified_names:
            simplified_names_list = simplified_names.split(", ")

        if "/" in simplified_names:
            simplified_names_list = simplified_names.split("/")

        if organism_name in simplified_names_list:
            target_species_value_matched_with_input = target_species_option.get_attribute("value")

    target_gene_element.send_keys(target_gene)
    target_species_element.select_by_value(target_species_value_matched_with_input)

    guide_rna_submit_button = driver.find_element_by_xpath("//*[@id=\"searchRequest\"]")
    # Submits form
    guide_rna_submit_button.click()

    try:
        wait = WebDriverWait(driver, time_allotment)
        wait.until(EC.presence_of_element_located((By.XPATH, "/html/body/main/div/div[3]/div/div/table")))

        guide_number = 1
        finished = False

        # Adds links of each guideRNA to obtain off target information
        off_target_site_links = []
        while not finished:
            try:
                guide_row_element = driver.find_element_by_xpath("//*[@id=\"guide" + str(guide_number) + "\"]")
                guide_rna_build_info = []

                for column in range(1, 12):
                    guide_row_and_column_cell_element = driver.find_element_by_xpath("/html/body/main/div/div[3]/div/div/table/tbody/tr[" + str(guide_number) + "]/td[" + str(column) + "]")
                    guide_rna_build_info.append(guide_row_and_column_cell_element.text)

                output.append(guide_rna_build_info)

                off_target_site_links.append(driver.current_url + "details/" + str(guide_number))

            except NoSuchElementException:
                finished = True

            guide_number += 1

        # Shifts through all guideRNA sites to obtain off target information and save it into
        # guide_rna_results_storage_folder
        guide_number = 1
        for off_target_site in off_target_site_links:
            driver.get(off_target_site)
            try:
                WebDriverWait(driver, time_allotment).until(EC.presence_of_element_located((By.XPATH, "/html/body/main/div/div[6]/div/div/table/tbody/tr[12]/td[1]")))
                table_elements = driver.find_elements_by_xpath("/html/body/main/div/div[5]/div/div/table/tbody/tr")

                # This section saves the off-target data into the guideRNA results folder
                # This save implementation is done differently than the other two, which employ
                # a simple CSV approach because the HI-CRISPR method required such folder output
                # instead of simple CSV output.
                num_rows = len(table_elements)
                off_target_info = ""
                off_target_row_info = []
                for table_element in table_elements:
                    td_elements = table_element.find_elements(By.TAG_NAME, "td")
                    if num_rows == 1 and len(td_elements) == 1:
                        off_target_info = "There are no predicted off-targets."
                    elif num_rows == 1 and len(td_elements) == 3:
                        off_target_info = ",%s,%s,%s" % (td_elements[0].text, td_elements[1].text, td_elements[2].text)
                    else:
                        off_target_row_info.append(",%s,%s,%s" % (td_elements[0].text, td_elements[1].text, td_elements[2].text))

                if off_target_info == "":
                    off_target_info = ";".join(off_target_row_info)

                guide_RNA_sequence = output[guide_number][1]
                guide_number_off_target_file = open("%s/%s.offtargets" % (off_target_folder, guide_number), "w")

                guide_number_off_target_file.write(guide_RNA_sequence + "\n")
                guide_number_off_target_file.write(off_target_info)
                guide_number_off_target_file.close()

                guide_number += 1

            except TimeoutException:
                driver.close()
                raise Exception("Finding off-target site information for knockout of " + str(target_gene) + " in organism " + str(organism_name) +" not found within " + str(time_allotment) + " seconds.")

    except TimeoutException:
        driver.close()
        raise Exception("guideRNAs for gene " + str(target_gene) + " in organism " + str(organism_name) + " not found within " + str(time_allotment) + " seconds.")

    driver.close()

    # Generates output summary of run
    guide_rna_output_summary_file = open("%s/Knockout Output Summary.txt" % (gene_folder), "w")
    guide_rna_output_summary_file.write("Date Run: %s" % (date.today()))
    guide_rna_output_summary_file.close()

    # Saves output, a 2D array, into csv format for future retrieval
    column_headers = output[0]
    non_column_headers = output[1:]

    with open("%s/guideRNAs.csv" % (gene_folder), "w") as f:
        write = csv.writer(f)

        write.writerow(column_headers)
        write.writerows(non_column_headers)

    return output

# Obtains archive guideRNA data and returns False if not found or corrupted
def obtainGuideRNAsAncientRunCHOPCHOP(organism, gene, guide_rna_results_storage_folder):
    if os.path.isdir("%s/%s/%s" % (guide_rna_results_storage_folder, organism, gene)):
        if os.path.isfile("%s/%s/%s/guideRNAs.csv" % (guide_rna_results_storage_folder, organism, gene)):
            df = pd.read_csv("%s/%s/%s/guideRNAs.csv" % (guide_rna_results_storage_folder, organism, gene), header = None)
            return df.values.tolist()
        else:
            return False
    else:
        return False

# Obtains age of guideRNA data and returns False if not found or corrupted
def obtainAgeOfSavedGuideRNAsinDays(organism, gene, guide_rna_results_storage_folder):
    if os.path.isdir("%s/%s/%s" % (guide_rna_results_storage_folder, organism, gene)):
        if os.path.isfile("%s/%s/%s/Knockout Output Summary.txt" % (guide_rna_results_storage_folder, organism, gene)):
            summary_info_file = open("%s/%s/%s/Knockout Output Summary.txt" % (guide_rna_results_storage_folder, organism, gene), 'r')
            raw_str_date = summary_info_file.read().replace('\n', '')
            summary_info_file.close()

            date_ran = raw_str_date.split(": ")[1]
            date_ran = datetime.strptime(date_ran, "%Y-%m-%d")

            date_today = datetime.now()

            return (date_today - date_ran).days
        else:
            return False
    else:
        return False

# Obtains guideRNA data if archived and younger than acceptable_age days, else fetches it from the CHOPCHOP website
def obtainGuideRNAsSavedOrOnlineCHOPCHOP(organism, gene, guide_rna_results_storage_folder, time_allotment = 120, selenium_driver = "chromedriver/chromedriver", acceptable_age = sys.maxsize):
    saved_guide_rnas = obtainGuideRNAsAncientRunCHOPCHOP(organism, gene, guide_rna_results_storage_folder)

    if not saved_guide_rnas:
        return obtainGuideRNAsOnlineCHOPCHOP(organism, gene, guide_rna_results_storage_folder, time_allotment = time_allotment, selenium_driver = selenium_driver)
    else:
        age_of_saved_run = obtainAgeOfSavedGuideRNAsinDays(organism, gene, guide_rna_results_storage_folder)

        if age_of_saved_run < acceptable_age:
            return obtainGuideRNAsAncientRunCHOPCHOP(organism, gene, guide_rna_results_storage_folder)
        else:
            return obtainGuideRNAsOnlineCHOPCHOP(organism, gene, guide_rna_results_storage_folder, time_allotment = time_allotment, selenium_driver = selenium_driver)
