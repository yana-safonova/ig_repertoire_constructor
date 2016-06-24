function launchSpec(dataProvider)
{
  var retval = {
    commandLine: ["python", "/ig_repertoire_constructor/BaseSpace/apps/runIGRC.py"],
    containerImageId : "docker.illumina.com/ig_repertoire_constructor/igrc:latest",
    Options: [ "bsfs.enabled=true" ]
  };

  return retval;
}
