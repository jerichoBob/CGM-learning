# gather
# this will serve as the tool to gather together all of the analysis distributed across the various folders within analysis/j1429.
#
# the first task is to gather together all of the PDF stackplots and put them into a single PDF file.

import sys, os
import argparse


def find_pdfs(basedir):
    return (os.path.join(root, file)
            for root, dirs, files in os.walk(basedir)
            for file in files if file.lower().endswith('.pdf'))


parser = argparse.ArgumentParser(
                    prog='pather.py',
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    description="""
This will serve as the tool to gather together all of the analysis distributed across the various folders within analysis/j1429.

Example usage:
    python gather.py -d analysis/j1429 -p
    >>> 

""", 
                    # epilog='----'
                    )
parser.add_argument('-d', '--basedir', required=True, help='the directory across which the gathering will occur')
parser.add_argument('-p', '--pdf',    action='store_true', help='gather together the PDF stackplots')

if __name__ == "__main__":

    args = parser.parse_args()

    if args.basedir is not None:
        basedir = args.basedir
        
        if args.pdf is not None:

            pdf_files = list(sorted(find_pdfs(basedir)))
            for pdf in pdf_files:
                print(pdf)
