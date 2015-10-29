"""
"""
import subprocess

output = subprocess.Popen(["python", "setup.py", "install", "--dry-run"],
                          stdout=subprocess.PIPE).communicate()[0]
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
print result
