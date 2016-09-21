bad_words=['OSI.insertDynamicalSystem']
import os

for dirname, dirnames, filenames in os.walk('.'):
    # print path to all subdirectories first.
    for subdirname in dirnames:
        print(os.path.join(dirname, subdirname))

    for filename in filenames[1:]:
        if ('.new' in filename):
            print(filename)
            os.rename(os.path.join(dirname, filename), os.path.join("./garbage/", filename))
        else:
            if ( '.py' in filename):
                print(os.path.join(dirname, filename))
                
                with open(os.path.join(dirname, filename)) as oldfile, open(os.path.join(dirname, filename+'.new'), 'w') as newfile:
                    for line in oldfile:
                        if not any(bad_word in line for bad_word in bad_words):
                            newfile.write(line)
                    
                os.rename(os.path.join(dirname, filename+'.new'), os.path.join(dirname, filename))

                
