# cellname

最近在弄的一个自动化注释的python，注意只能做人和小鼠的
其他物种需要转换symbol，同时适配空间组和单细胞，反正输入的是anadata格式，
什么类型的数据并不care。

不过我自己测试的感觉，单细胞用`q = 0.5`,空间组用`q = 0.75`比较好

安装就是

```angular2html
git clone https://github.com/wangjiaxuan666/cellname
cd cellname
pip install -e .
```

## Tutorial

见[NOTEBOOK/Tutorials.ipynb](notebook/Tutorials.ipynb)
