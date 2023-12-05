# gather
# this will serve as the tool to gather together all of the analysis distributed across the various folders within analysis/j1429.
#
# the first task is to gather together all of the PDF stackplots and put them into a single PDF file.

import os, re
import argparse
from PyPDF2 import PdfWriter, PdfReader, PageObject, PaperSize, Transformation
from PyPDF2.generic import AnnotationBuilder

from pathlib import Path
from typing import Union, Literal, List
from datetime import datetime

from bobutils import utils as bu


def find_pdfs(basedir):
    return (os.path.join(root, file)
        for root, dirs, files in os.walk(basedir)
            for file in files if file.lower().endswith('.pdf'))

def print_list(list):
    for pdf in list:
        print(pdf)

def paths_to_tree(paths):
    """Convert a list of paths to a nested dictionary tree."""
    tree = {}
    for path in paths:
        parts = path.split(os.sep)
        node = tree
        for part in parts:
            node = node.setdefault(part, {})
    return tree

def print_tree(tree, indent='', prefix=''):
    """Pretty print the tree."""
    items = list(tree.items())
    for i, (key, value) in enumerate(items):
        is_last = i == len(items) - 1
        if value:
            print(f"{indent}{prefix}── {key}")
            next_prefix = '└─' if is_last else '├─'
            next_indent = indent + ('    ' if is_last else '│   ')
            print_tree(value, indent=next_indent, prefix=next_prefix)
        else:
            print(f"{indent}{prefix}── {key}")


# output_file_page_number = 0

def label_pages(label: str, content_pdf: Path, writer: PdfWriter):
    
    reader = PdfReader(content_pdf)
    page_indices = list(range(0, len(reader.pages)))

    # Create the annotation and add it
    xLL = 10
    yLL = 10
    width = 600
    height = 20
    annotation1 = AnnotationBuilder.free_text(label,
                                             rect=(xLL, yLL, xLL+width, yLL+height), # [xLL, yLL, xUR, yUR]
                                             font="Arial",
                                             bold=True,
                                             italic=True,
                                             font_size="6000pt",
                                             font_color="ff0000",
                                             border_color="000000",
                                            #  background_color="cdcdcd",
    )
    base_page_number = len(writer.pages)
    user_units = 72.0 # 72 points per inch
    target_page_width = float(11.0 * user_units)
    target_page_height = float(8.5 * user_units)

    for index in page_indices:

        content_page = reader.pages[index]     
        orig_page_width = float(content_page.artbox.width)
        orig_page_height = float(content_page.artbox.height)

        scale_factor_w = target_page_width / orig_page_width
        scale_factor_h = target_page_height / orig_page_height
        # print(f"scale factor w: {bu.sig_figs(scale_factor_w,5)}  scale factor h: {bu.sig_figs(scale_factor_h,5)}")
        scale_factor = min(scale_factor_w, scale_factor_h) / 1.0
        scaled_page_width = orig_page_width * scale_factor
        scaled_page_height = orig_page_height * scale_factor
        # print(f"scale factor: {bu.sig_figs(scale_factor,5)}")
        tx = (target_page_width - scaled_page_width) / 2.0
        ty = (target_page_height - scaled_page_height) / 2.0

        # Scale the page
        content_page.add_transformation(Transformation().scale(scale_factor))

        # print(f"orig page width:   {bu.sig_figs(orig_page_width/72.0, 5)  }  orig page height:   {bu.sig_figs(orig_page_height/72.0, 5)}")
        # print(f"scaled page width: {bu.sig_figs(scaled_page_width/72.0, 5)}  scaled page height: {bu.sig_figs(scaled_page_height/72.0,5)}")

        new_page = PageObject.create_blank_page(None, target_page_width, target_page_height)
        new_page.merge_page(content_page)
        new_page.add_transformation(Transformation().translate(tx=tx, ty=ty))

        writer.add_page(new_page)
        writer.add_annotation(page_number=base_page_number+index, annotation=annotation1)

def get_sightline_from_path(path: str) -> Union[str, None]:
    """Extract the sightline from the path."""
    sightline_match = re.search(r'[\w\d\_\+]+/(\d+)_', path)
    if sightline_match:
        sightline = sightline_match.group(1)
    else:
        sightline = None
    return sightline
def get_redshift_from_path(path: str) -> Union[float, None]:
    """Extract the redshift from the path."""
    redshift_match = re.search(r'z_([\d.]+)', path)
    if redshift_match:
        redshift = float(redshift_match.group(1))
    else:
        redshift = None
    return redshift
def get_ion_index_from_path(path: str) -> Union[int, None]:
    """Extract the ion index from the path."""
    ion_index_match = re.search(r'Ions(\d+)', path)
    if ion_index_match:
        ion_index = int(ion_index_match.group(1))
    else:
        ion_index = None
    return ion_index

def gather_pdfs(paths):
    """Gather all PDFs into a single PDF file. Label each page with the sightline and redshift"""
    writer = PdfWriter()

    # now = datetime.utcnow()
    # timestamp = now.isoformat() + "Z"
    # output_file = "./output-" + timestamp + ".pdf"

    output_file = "./output.pdf"
    # output_file_page_number = 0

    with open(output_file, "wb") as fp:
        for filepath in paths:
            label = f"J1429 - Sightline: #{get_sightline_from_path(filepath)}  z: {get_redshift_from_path(filepath)}  Ion Grouping: {get_ion_index_from_path(filepath)}"
            print("- "*20)
            # print(f"filepath: {filepath}")
            print(f"   label: {label}")
            label_pages(label, filepath, writer)
            writer.write(fp)

    writer.close()



parser = argparse.ArgumentParser(
                    prog='pather.py',
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    description="""
This tool gathers together the analysis distributed across various folders within `basedir``.

Example usage:
# gather all PDFs and pretty print the list of gathered PDFs
    python gather.py -d analysis/j1429 -p -pp
    
# gather all PDFs and merge them into a single PDF file     
    python gather.py -d analysis/j1429 -p -m 

""", 
                    # epilog='----'
                    )
parser.add_argument('-d', '--basedir', required=True, help='the directory across which the gathering will occur')
parser.add_argument('-p', '--pdf',    action='store_true', help='gather together the PDF stackplots')
parser.add_argument('-m', '--merge',    action='store_true', help='merge gathered PDF stackplots into a single PDF file')
parser.add_argument('-pp', '--prettyprint',  action='store_true', help='print out the list of gathered PDF stackplots in a pretty format')

if __name__ == "__main__":

    args = parser.parse_args()
    print(f"args: {args}")
    if args.basedir is not None:
        basedir = args.basedir

        if args.pdf:

            pdf_files = list(sorted(find_pdfs(basedir)))

            if args.prettyprint:
                tree = paths_to_tree(pdf_files)
                print_tree(tree)

            if args.merge:
                print("merging PDFs...")
                # print_list(pdf_files)
                gather_pdfs(pdf_files)
            
            
