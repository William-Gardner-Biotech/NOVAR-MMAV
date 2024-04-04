# OVerall aVg coverAge foLd

import logging
import argparse
from datetime import datetime

# Create a unique filename with the current date and time for debug
current_time = datetime.now().strftime("%m-%d_%H-%M-%S")
filename = f"debug_{current_time}.log"

logging.basicConfig(filename=filename, level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

def get_args():
    parser = argparse.ArgumentParser(description='Calculate Overall aVg coverAge foLd value for each mpileup file')

    # Input file argument
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='The path to the input mpileup file')

    arguments = parser.parse_args()

    return arguments

def get_oval(mpileup_file):
    """
    To get Overall aVg coverAge foLd for this we must take the avg_fold,
    multiply it by the % genome coverage to get the  avg coverage depth per segment.
    Then take that number and relate it with the segment length. Combine those numbers
    and rebuild it to be the whole genome % genome coverage depth.
    """
    with open(mpileup_file, 'r') as cov_file:
        cov_file.readline()
        whole_genome_length = 0
        deconvoluted = 0
        for line in cov_file:
            values = line.split('\t')
            # values[0]:id, [1]:avg_fold, [2]:length, [4]:covered_percent
            val = float(values[1])*float(values[4])/100
            logging.debug(f"{values[0].split(' ')[0]} aVg coverAge foLd {val}")
            whole_genome_length += int(values[2])
            # This will give us the segment reads
            deconvoluted += val*int(values[2])
        
        oval = deconvoluted/whole_genome_length
        logging.debug(f"Overall aVg coverAge foLd {oval:.2f}X")

        # Print statement will send the value to the terminal which we can pickup in linux
        print(f"{oval:2f}")
        

def main():
    args = get_args()
    get_oval(args.input)

main()