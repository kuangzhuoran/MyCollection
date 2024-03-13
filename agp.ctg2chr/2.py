#!/usr/bin/python3


#如果有一个'目录',内容是          A  
#                                 C
#                                 D
#                                
#同时有一个'库'，其内容包含了'目录'，
#具体为                 XXX A YYY
#                       XXX B ZZZ
#                       XXX C PPP
#                       XXX D PPP
#                       XXX E PPP
#这种格式
#
#把'库'能匹配到'目录'的行筛选出来
#换言之，'库'的内容我只取A、C、D这三行
#就可以用这个脚本
#最后得到 XXX A YYY
#         XXX C PPP
#         XXX D PPP
#


import sys
import re

list1 = {}
#f = open("tmp.sam",'r')  #待筛选的文件-库
f = open(sys.argv[1],'r') #文件多的话就用传入参数
for line in f:
    line = line.strip()
    content = line.split('\t')  #具体分隔符，具体第几列
    name = content[0]
    list1[name] = line
    
#freqfile = "tmp.ctg"  #要从库里面筛选出的东西-目录
#f = open(freqfile, 'r')  #这里没用传入参数

f = open(sys.argv[2],'r') #用传入参数就这样
for line in f:
    line = line.strip()
    content1 = line.split('\t')
    name1 = content1[0]

    if name1 in list1:
        print(line + "\t" + list1[name1])
        #print(list1[name1])
