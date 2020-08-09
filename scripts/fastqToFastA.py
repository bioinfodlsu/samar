fp = open(snakemake.input[0])
fpW = open(snakemake.output[0],"a")
count = 1
readCount = 1
for i, line in enumerate(fp):
    if(i == count):
    	fpW.write(">_Read_%d/%s\n" % (readCount,snakemake.output[0].split("_")[1].split(".")[0]))
    	fpW.write(line)
    	readCount+=1
    	count+=4
fp.close()
fpW.close()

