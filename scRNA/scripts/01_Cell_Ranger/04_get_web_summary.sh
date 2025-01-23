#!/bin/bash

# source environment variables
source ../../refs/.env

# extract webs summaries
for sample in E3CF1 E3CF2 E3CM1 E3CM2 E3PF1 E3PF2 E3PM1 E3PM2 E4CF1 E4CF2 E4CM1 E4CM2 E4PF1 E4PF2 E4PM1 E4PM2
do
	cd $PROJECT_DIR/counts/$sample/outs
	cp web_summary.html ../../web_summaries/"$sample"_web_summary.html
done
