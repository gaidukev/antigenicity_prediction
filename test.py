from pickle import load, dump

noise_file = open("noise", "rb")
record = load(noise_file)
noise_file.close()
'''
for k,v in record.items():
    if v >= 4 and k[-1] == "A":
        print(k+":", v)
'''
import winsound

winsound.Beep(2500, 1000)
