import os

def getLatest(path) :
  latest = -1
  for filename in os.listdir(path) :
    filenames = filename.split('-')
    if len(filenames) > 2 :
      try:
        i = int(filenames[0])
      except ValueError:
        i = 0
    if latest < i :
      latest = i
  return str(latest)

linux_intel_latest = getLatest('/home/autotest/buildmaster/intel/OpenCMISS')
linux_gnu_latest = getLatest('/home/autotest/buildmaster/gnu/OpenCMISS')

text = 'ALIASES += linux_intel_latest="'+linux_intel_latest+'"\nALIASES += linux_gnu_latest="'+linux_gnu_latest+'"\n'

file =open('Doxyfile', 'a')
file.write(text)
file.close()

os.system("doxygen Doxyfile")

