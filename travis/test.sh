
#!/bin/bash
cromwell_version="33"
wdl="processing-for-variant-discovery-gatk4.wdl"
json="processing-for-variant-discovery-gatk4.json"

#downloand cromwell
wget http://github.com/broadinstitute/cromwell/releases/download/$cromwell_version/cromwell-$cromwell_version.jar
#excute wdl script
java -Dconfig.file=travis/resources/google-adc.conf -Dbackend.providers.JES.config.project=broad-dsde-outreach -Dbackend.providers.JES.config.root=gs://beri-test/gatk-workflows-travis-tests/ -jar cromwell-$cromwell_version.jar run $wdl -i $json
