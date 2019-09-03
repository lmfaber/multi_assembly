import logging
from os.path import basename
from os.path import splitext

inputFiles = snakemake.input
outputFiles = snakemake.output

# Create logger
FORMAT = '{asctime} [{levelname}]: {message}'
logging.basicConfig(level=logging.INFO, format=FORMAT, style='{', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

for inputFile, outputFile in zip(inputFiles, outputFiles):
    baseName = splitext(basename(inputFile))[0]
    logger.info(f'Renaming header for contigs from {baseName} and saving to {outputFile}')
    with open(inputFile, 'r') as reader, open(outputFile, 'w') as writer:
        i = 1
        for line in reader:
            if line.startswith('>'):
                baseName = baseName.upper()
                line = f'>{i}_{baseName}\n'
                i += 1
            writer.write(line)