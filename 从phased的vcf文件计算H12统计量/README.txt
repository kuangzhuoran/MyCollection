#Usage:
#Step1
python  1.calc_garud_h_from_vcf.py  -i test.vcf  -w 20000  -s 5000  --min-sites 1  -o test.20k_5k.output.txt
python  1.calc_garud_h_from_vcf_scikit_allel.py  -i test.vcf  -w 20000  -s 5000  --min-sites 1  -o test.20k_5k.output.txt

#Step2
python 2.normalize_h2h1_by_h12.py -i test.20k_5k.output.txt -o test.20k_5k.output.normalized.txt --clip

#这里的脚本或者说明都是ChatGPT生成的
#请详细阅读脚本说明或者脚本头的注释

#当前版本是v20260401
