import logging

inputFiles = snakemake.input
if len(snakemake.output) == 2:
    outputFiles = [snakemake.output['R1'], snakemake.output['R2']]
else:
    outputFiles = [snakemake.output['R1']]

# Create logger
FORMAT = '{asctime} [{levelname}]: {message}'
logging.basicConfig(level=logging.INFO, format=FORMAT, style='{', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

def createNewHeader(line):
    line = line.split('/')
    firstPart = line[0]
    firstPart = firstPart.replace('_R', '')
    firstPart = firstPart.replace('_F', '')
    memberOfPair = line[-1].replace('\n', '')
    newLastPart = '{}:N:0:ACACAC'.format(memberOfPair)
    return('{} {}\n'.format(firstPart, newLastPart))

for inFile, outFile in zip(inputFiles, outputFiles):
    logger.info('Processing: {}'.format(inFile))
    try:
        with open(inFile, 'r') as reader, open(outFile, 'w') as writer:
            row = 1
            for line in reader:
                if row == 1:
                    newHeader = createNewHeader(line)
                    writer.write(newHeader)
                elif row == 3:
                    writer.write('+\n')
                else:
                    writer.write(line)
                
                
                if row == 4:
                    row = 1
                else:
                    row += 1
    except:
        logger.error('Failed the processing of this file: {}'.format(inFile))