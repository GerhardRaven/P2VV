#pdftk A="${INFILE}" cat A"${PAGE}" output "${OUTFILE}"
#pdfcrop --margins "3 3 3 3" "${OUTFILE}"
## Times-Roman --> LinBiolinum
#sed -i -e 's/Times-Roman/LinBiolinum/g' "${OUTFILE}"
#gs -sFONTPATH=/project/bfys/jleerdam/local/root/fonts -o "${OUTFILE}" -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress "${OUTFILE}"
#gs -o "${OUTFILE}" -sDEVICE=pdfwrite -sProcessColorModel=DeviceCMYK -sColorConversionStrategy=CMYK -sColorConversionStrategyForImages=CMYK "${OUTFILE}"

import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '--inputFile',  '-i', required = True )
parser.add_argument( '--outputFile', '-o' )
parser.add_argument( '--pages', '-p', nargs = '+', type = int, default = [ 1 ] )
parser.add_argument( '--fontPath', '-f', default = 'fonts' )
parser.add_argument( '--fontName', '-n', default = 'LinBiolinum' )
parser.add_argument( '--margins', '-m', nargs = '+', type = int, default = [ 3, 3, 3, 3 ] )

args = parser.parse_args()
assert args.inputFile.endswith('.pdf')
outFile = args.outputFile if args.outputFile else args.inputFile
assert outFile.endswith('.pdf')
assert len( args.pages ) > 0
assert len( args.margins ) == 4

import os, subprocess, tempfile
tmp = tempfile.NamedTemporaryFile()
tmpFile = os.path.realpath(tmp.name)
tmp.close()
for page in args.pages :
    subprocess.call( [ 'pdftk', 'A=%s' % args.inputFile, 'cat', 'A%d' % page, 'output', tmpFile ] )
    subprocess.call( [ 'pdfcrop', '--margins', '%d %d %d %d' % tuple( args.margins ), tmpFile, tmpFile ] )
    subprocess.call( [ 'sed', '-i', '-e', 's/Times-Roman/%s/g' % args.fontName, tmpFile ] )
    subprocess.call( [ 'gs', '-sFONTPATH=%s' % args.fontPath, '-sDEVICE=pdfwrite', '-dPDFSETTINGS=/prepress', '-o', outFile
                      , tmpFile ] )
    subprocess.call( [ 'rm', tmpFile ] )
