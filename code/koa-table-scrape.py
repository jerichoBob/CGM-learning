# imports
import requests
from bs4 import BeautifulSoup
import sys, os
import argparse


parser = argparse.ArgumentParser(
                    prog='koa-table-scrape.py',
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    description="""
Given a URL to a KOA table, this script will create a CSV file of the table data.
https://koa.ipac.caltech.edu/cgi-bin/bgServices/nph-bgExec?bgApp=/KOA/nph-KOA&instrument_kc=kcwi&filetype=science&calibassoc=assoc&locstr=14%3A29%3A54+%2B12%3A02%3A35&regSize=30&resolver=ned&radunits=arcsec&spt_obj=spatial&single_multiple=single
""")
parser.add_argument('--url', required=True, help='the URL to the page that holds the KOA table')

if __name__ == "__main__":
    args = parser.parse_args()
    if args.url is None:
        parser.print_help(sys.stderr)
        sys.exit(1)
    else: 
        URL = args.url
        page = requests.get(URL)
        soup = BeautifulSoup(page.content, "html.parser")
        print(soup.prettify())
        # table = soup.find(id="/html/body/table[2]/tbody/tr/td/div")
        # print(table.prettify())



    # main()