function launchSpec(dataProvider)
{
  var retval = {
    commandLine: ["python", "/opt/ig_repertoire_constructor/src/extra/BaseSpace/apps/runIGREC.py"],
    containerImageId : "docker.illumina.com/ig_repertoire_constructor/igrec:latest",
    Options: [ "bsfs.enabled=true" ]
  };

  return retval;
}
