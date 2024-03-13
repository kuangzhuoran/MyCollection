#!bash
file=$1
window=$2
help() {
    sed -rn 's/^### ?//;T;p;' "$0"
}

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
            help
            exit 1
fi

####file:第二列是位置,第三列是要累加的数值。
####bash window.sh  file    10000
####只支持单条染色体来，多条请用window2.sh
####数据第二列最大值决定多少个窗的输出，比如染色体100M，file最后只有80M，那么造窗直到80M

sort -k2,2n  --parallel=10  $file | 
awk '
        BEGIN{
                sum=0;  #窗里面多少个点
                ehh=0;
        };

        {       
                pos[$2]=$2;  ##  $2 pos所在的列
		maxp[NR]=$2; ###记录最大的位点
                xpehh[$2]=$3;  ##  $3  统计的东西的列
        };
	END{	
		fornum=(maxp[NR]/'$window')+1
		for (i=0;i<=fornum;i++){   
			l=int(i*'$window');
			r=int((i+1)*'$window'); 
			for(a in pos){    
				p=int(a)
				if (p>=l && p<=r){
                  	              sum++; ##有多少个点
                                      ehh+=xpehh[p]; ###统计的东西 以ehh 为例
	                        }
			}
			if(sum>0){
				print l"\t"r"\t"ehh"\t"sum;
			}
			if(sum==0){
				print  l"\t"r"\t0\t0"
			}
			sum=0;
			ehh=0;
		}
        }
' 

