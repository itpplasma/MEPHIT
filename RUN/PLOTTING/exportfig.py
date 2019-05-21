import matplotlib.pyplot as plt
import os

def exportfig(filename):
    print('export: ' + filename)
    plt.tight_layout()
    plt.savefig(filename + '.pdf')
    os.system('pdfcrop ' + filename + '.pdf ' + filename + '.pdf')
    
