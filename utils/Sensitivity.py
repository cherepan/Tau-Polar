#!/usr/bin/env python

import os
import argparse
import ROOT
import math
import array



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--input-file",help="input file. [Default: %(default)s] ", action="store", default = 'TauolaHelicity_results.root')
    parser.add_argument("-r", "--obs",help="Observable [Default: %(default)s] ",  type=str, action="store", default = 'omega_a1')
    parser.add_argument('-p', '--pol-value',help="Polarization value [Default: %(default)s] ", default=0, )
    args = parser.parse_args()



    def get_histo(path):
        input_file_path, hist_name = path.split(":")
        input_file = ROOT.TFile(input_file_path, "READ")
        hist = input_file.Get(hist_name)
        hist.SetDirectory(0)
        input_file.Close()
        return hist

    histoplus = get_histo(args.input_file+":"+args.obs+"_plus");
    histominus = get_histo(args.input_file+":"+args.obs+"_minus");
    

    histoplus.Scale(1./histoplus.Integral())
    histominus.Scale(1./histominus.Integral())

    hf = histoplus.Clone("hf")
    hg = histoplus.Clone("hg")
    hf.Add(histominus)
    hg.Add(histominus,-1)

    hf.Scale(0.5)
    hg.Scale(0.5)

    hf.Add(hg,args.pol_value)
    hg.Multiply(hg)

    hg.Divide(hf)


#    histoplus.Draw()
#    histominus.Draw("same")
    print "sensitivity =", math.sqrt(hg.Integral())
