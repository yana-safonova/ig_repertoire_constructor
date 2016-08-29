from sample_reader import SampleReader

samples_fv = SampleReader().read('AbVitro/flu_time_course/FV/', ['25'])
samples_gmc = SampleReader().read('AbVitro/flu_time_course/GMC/', ['8'])
samples_ido = SampleReader().read('AbVitro/flu_time_course/IDO/')
samples_age = SampleReader().read('age/')
samples_paired = SampleReader().read('AbVitro/paired/')
