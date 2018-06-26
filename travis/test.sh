
#!/bin/bash
#declare wdls to execute

#downloand cromwell
wget http://github.com/broadinstitute/cromwell/releases/download/32/cromwell-32.jar
#excute java script
java -Dconfig.file=travis/resources/google-adc.conf -Dbackend.providers.JES.config.project=broad-dsde-outreach -Dbackend.providers.JES.config.root=gs://beri-test/ -jar cromwell-32.jar run processing-for-variant-discovery-gatk4.wdl -i processing-for-variant-discovery-gatk4.hg38.wgs.inputs.json
