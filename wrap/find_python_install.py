"""
"""
import subprocess
import locale

encoding = locale.getdefaultlocale()[1]

output = subprocess.Popen(["@Python3_EXECUTABLE@", "setup.py", "--dry-run", "install"],
                          stdout=subprocess.PIPE).communicate()[0]

if encoding:
    output = output.decode(encoding).split('\n')
else:  # let's cross fingers here ...
    output = output.split('\n')

for line in output:
    if line.count('egg') and line.count('fake'):
        result = line
        break
result = result.split(' ')
for line in result:
    if line.count('egg'):
        result = line
        break
result = result.split('fake')[0]
print(result)
