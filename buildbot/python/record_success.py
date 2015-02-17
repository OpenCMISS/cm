import os, subprocess
from jinja2 import Template, Environment, FileSystemLoader

opencmiss_repository = os.environ['OPENCMISS_ROOT']

os.chdir("%s/cellml"%(opencmiss_repository))
proc_cellml = subprocess.Popen(["git rev-list HEAD --max-count=1"], stdout=subprocess.PIPE, shell=True)
(out_cellml, err_cellml) = proc_cellml.communicate()

os.chdir("%s/cm"%(opencmiss_repository))
proc_cm = subprocess.Popen(["git rev-list HEAD --max-count=1"], stdout=subprocess.PIPE, shell=True)
(out_cm, err_cm) = proc_cm.communicate()

os.chdir("%s/examples"%(opencmiss_repository))
proc_examples = subprocess.Popen(["git rev-list HEAD --max-count=1"], stdout=subprocess.PIPE, shell=True)
(out_examples, err_examples) = proc_examples.communicate()

f = open("%s/build/logs/last_successful_build.html"%(opencmiss_repository), "w")
# configure jinja2 environment
env = Environment(loader=FileSystemLoader('.'), extensions=['jinja2.ext.loopcontrols'])

# main code
os.chdir("%s/cm/buildbot"%(opencmiss_repository))
t = env.get_template("template/last_success_build.html")
f.write(t.render(cellml_commits=out_cellml,
                 cm_commits=out_cm,
                 examples_commits=out_examples))


