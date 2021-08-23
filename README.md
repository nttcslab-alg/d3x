# d3x

This is an implementation of algorithm d3x for finding all the solutions of [exact cover problems](https://en.wikipedia.org/wiki/Exact_cover). d3x accelerates [Algorithm DLX](https://en.wikipedia.org/wiki/Dancing_Links) by comporessing the input by using [Zero-suppressed Binary Decision Diagram](https://en.wikipedia.org/wiki/Zero-suppressed_decision_diagram). For details, please refer to the paper.

## requirements
- c++ compiler supporting c++17 (gcc, clang)
- cmake 3.16
## compile

```bash
$ mkdir build
$ cd build
$ cmake ..
$ cmake -build .
```

## run

```bash
$ ./d3x -z zdd_file
```
- `zdd_file` follows the format of ZDDs used in [Graphillion](https://github.com/takemaru/graphillion).  In graphillion, you can get the ZDD corresponding to a GrpahSet object by using `gs.dump(fp)` method.

## Reference

Masaaki Nishino, Norihito Yasuda, and Kengo Nakamura, "Compressing Exact Cover Problems with Zero-suppressed Binary Decision Diagrams", in Proc. of the 30th International Joint Converence on Artificial Intelligence (IJCAI 21), [Paper](https://www.ijcai.org/proceedings/2021/275) 



