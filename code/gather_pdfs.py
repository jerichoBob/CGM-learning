from PyPDF2 import PdfWriter, PdfMerger, PdfFileReader, PdfReader, PageObject, PaperSize, Transformation
import os, re, argparse

def find_pdfs(basedir):
    return (os.path.join(root, file)
        for root, dirs, files in os.walk(basedir)
            for file in files if file.lower().endswith('.pdf'))
                # if re.match(r"Lecture\d+.pdf", file))
                # if file.lower().endswith('.pdf'))

def print_list(list):
    for pdf in list:
        print(pdf)

def merge_pdfs(basedir, outfile):
    """Gather all PDFs from {basedir} into a single PDF file {outfile}. A simple copy/paste exercise"""

    pdf_files = list(sorted(find_pdfs(basedir)))
    print(f"pdf_files in {basedir}")
    print_list(pdf_files)

    fh_output = PdfMerger()

    for filename in pdf_files:
        fh_input = open(filename, 'rb')
        reader = PdfReader(fh_input)
        fh_output.append(reader)
    
    fh_input.close()
    fh_output.write(outfile)


app_name = "gather_pdfs.py"
parser = argparse.ArgumentParser(
                    prog=app_name,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    description=f"""
{app_name} merges a sorted list of all pdfs contained in --basedir into a single output.pdf.

Example usage:
    python ~/bin/gather_pdfs.py -d . -o "Numerical_Methods_in_Physics_with_Python.pdf" 
    >>> outputs a file called output.pdf

""", 
                    # epilog='----'
                    )
parser.add_argument('-d', '--basedir', required=True, help='the source directory of the pdf files')
parser.add_argument('-o', '--outfile', help='override default output file name (default: output.pdf)')

if __name__ == "__main__":

    args = parser.parse_args()
    print(f"args: {args}")
    if args.basedir is not None: basedir = args.basedir
    else: 
        raise ValueError('gather_pdfs: basedir is None! Check your command line arguments')
    if args.outfile is not None: outfile = args.outfile
    else: outfile = "output.pdf"
    
    print(f"basedir: {basedir} outfile: {outfile}")
    merge_pdfs(basedir, outfile)
            
            
