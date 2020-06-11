<!--
#markdown source of Mapoet Niphy 
#
#Author    :  Mapoet Niphy
#Date      :  2019
#Institude :  SHAO
#
-->
# GEOTOOL

在测绘、地球科学、空间科学、大气科学以及近地行星科学研究，时间与空间系统相互转换是一切研究的基础，为此提供了相关的时间系统/时间格式，空间系统/空间格式的源码，及通过该源码实现的时间或者坐标转化小工具。该小工具支持以流的形式将时间或者坐标序列转化为需求的形式。

## 功能简介

著名的[SOFA软件](http://www.iausofa.org/)，提供了天球坐标系统与地固坐标系相互转函数，也提供了不同时间系统的时间转化函数。我们利用了SOFA中的MJD与YMD转化模块，实现了以下时间形式转化：
* (1) YYYY/MM/DD [H:M:S],
* (2) YYYY/DOY [H:M:S],
* (3) MJD,
* (4) WEEK/DAY [SOD(Sec of day)] in GPS,
* (5) WEEK/SOW(Sec of Week) in GPS。

坐标形式主要实现了：
* (1) 大地坐标,
* (2) 极坐标,
* (3) 笛卡尔坐标。

而坐标系统基于[XFORM](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/sxform.html)，支持地球惯性坐标系（gei），地理坐标系（geo），地心太阳坐标系（gse），地心太阳磁坐标系（gsm），地磁坐标系（mag），太阳磁坐标系（sm）间的转换。

## 使用方式

以下介绍两个工具`DATETIMES`及`POSCONVERT`的使用方式。

### `DATETIMES`

在`Linux`可以通过，`echo`及`cat`命令将数据传入到执行程序中，而在`Windows`需要指向输入符来进行处理，基本调用形式如下：

* `Windows`
`DATETIMES abc < inputs > outputs`
* `Linux`
`cat inputs |./DATETIMES abc > outputs`
或`./DATETIMES abc <inputs > outputs`
一般可以直接使用`echo dates|./DATETIMES abc`，以嵌入`bash`命令中灵活使用。

以上`abc`指代转换参数, 可以是$$t_i t_j s$$或$$t_i o s$$，前者是形式转化，后者是时间计算。
$$t_i, t_j$$为时间格式标签，见时间形式及其前面的编号。
$$o$$表示计算方式，支持:(0),加；(9)，减；(8),时间差。
$$s$$为时间是否含有`time`部分，有以下含义：(0),只有日期；(1)，含有时分秒。

### `POSCONVERT`

类似于`DATETIMES`，`POSCONVERT`在`Windows`或`Linux`使用也有些许区别，具体如下：
* `Windows`
`POSCONVERT a[b]　[`$$s_i$$`2`$$s_j$$`:time] < inputs > outputs`
* `Linux`
`cat inputs |./POSCONVERT a[b]  [`$$s_i$$`2`$$s_j$$`:time] > outputs`
或`./POSCONVERT a[b] [`$$s_i$$`2`$$s_j$$`:time] <inputs > outputs`
一般可以直接使用`echo pos|./DATETIMES a[b] [`$$s_i$$`2`$$s_j$$`:time]`，以嵌入`bash`命令中灵活使用。

以上`a[b]`指代转换形式，可以是上面已列的坐标形式；`b`表示不同于`a`的形式，缺省`b`的情况表示坐标形式不变。` [`$$s_i$$`2`$$s_j$$`:time] `表示转入转出的坐标系统与参考时间，缺省情况下表示只对坐标形式转换。

由于坐标与时间密切相关，特别是惯性坐标系与地固坐标系，或者是ｘ轴指向为春分点或是地日连线的区别，以及地磁两极位置的变化等这些对坐标的影响比较大，所以在进行不同坐标系统转换时必须提供时间`time`；这里的坐标系支持所列的地球惯性坐标系（gei），地理坐标系（geo），地心太阳坐标系（gse），地心太阳磁坐标系（gsm），地磁坐标系（mag），太阳磁坐标系（sm），而时间格式支持包含时分秒的时间形式。

## 注意事项及使用案例　

注意事项，在坐标转换过程中，大地坐标形式只适合于地理坐标系，而其他两种在所有坐标系统都适合。但是程序中并没有作此判断，使用者需注意这方面问题。
具体使用案例见[时间及坐标转换样例](https://github.com/Mapoet/GEOTOOL/blob/master/TEST/T_DATETIMES%2BPOS.sh)。