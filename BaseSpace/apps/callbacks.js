function launchSpec(dataProvider)
{
  var retval = {
    commandLine: ["python", "/ig_repertoire_constructor/BaseSpace/apps/runIGRC.py"],
    containerImageId : "eodus/mydock",
    Options: [ "bsfs.enabled=true" ]
  };

  return retval;
}
