# CellName

最早是发现了一个比较好的Cell Marker 数据库，后面发现这个数据库的作者还开发一个自动化注释的软件
叫[adobo](https://github.com/oscar-franzen/adobo/tree/master/adobo)

因为我又测试过数据库的marker，感觉还挺准的，所以就想着这个注释的软件应该也不错。
但其实我发现这个软件问题还挺多的。

- 三年没更新了，放弃维护了，所以有个BUG根本跑不通
- adobo想做成一个和scanpy一样的单细胞分析框架，所以有不少功能可能三年前很OK，现在完全没必要了
- 其次作者没开发完全，有些地方逻辑不好，并只能做小鼠的
- 把参数写死了，没有灵活性

SO，自己动手丰衣足食，熬夜修BUG，平常工作间隙，偶尔抽个时间，写几下，完工！

目前只能做人和小鼠的，其他物种需要转换symbol，同时适配空间组和单细胞，反正输入的是anadata格式，
什么类型的数据并不care。

另外，我自己测试的感觉，单细胞用`q = 0.5`,空间组用`q = 0.75`比较好

安装就是

```angular2html
git clone https://github.com/wangjiaxuan666/cellname
cd cellname
pip install -e .
```

## 附带个教程，给项目组的师弟师妹们

## Tutorial

见[NOTEBOOK/Tutorials.ipynb](notebook/Tutorials.ipynb)
