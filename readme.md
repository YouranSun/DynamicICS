## Efficient Influential Community Discovery over Dynamic Graphs

### Compile

```bash
cmake -S ./
make -j4
```

### Run

#### Build

```bash
./bin/ICSDynamic <graph> <index> <m2> b advanced
```

* `<graph>`: input path - a undirected graph dataset, in the following format:

  ```
  n m
  u_1 v_1
  ...
  u_m v_m
  w_1 ... w_m
  ```
* `<index>`: output path - the ICD-Order index.
* `<m2>`: the index is built with `m - m2` among all the `m` edges in the graph.

#### Insertion / `OrdIns`

```bash
./bin/ICSDynamic <index-in> <index-out> <m2> i advanced <logfile>
```

* `<index-in>`: input path - the ICD-Order index to update.
* `<index-out>`: output path - the ICD-Order index after updating.
* `<m2>`: the number of edges to insert.
* `<logfile>`: log file path.

#### Deletion / `OrdDel`

```bash
./bin/ICSDynamic <index-in> <index-out> <m2> d advanced <logfile>
```

#### Query

```bash
./bin/ICSDynamic <index> <result> 0 q advanced <query>
```

* `<index>`: input path - the index to query
* `<result>`: output path - the top-$r$ $k$-ICs.
* `<query>`: input path - the queries, in the following format:

  ```
  q
  k_1 r_1
  ...
  k_q r_q
  ```
