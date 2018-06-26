
#!/bin/bash

#downloand cromwell
wget http://github.com/broadinstitute/cromwell/releases/download/33/cromwell-33.jar
#excute wdl script
java -Dconfig.file=travis/resources/google-adc.conf -Dbackend.providers.JES.config.project=broad-dsde-outreach -Dbackend.providers.JES.config.root=gs://beri-test/gatk-workflows-travis-tests/ -jar cromwell-32.jar run processing-for-variant-discovery-gatk4.wdl -i processing-for-variant-discovery-gatk4.hg38.wgs.inputs.json
