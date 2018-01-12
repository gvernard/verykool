#ifndef TABLE_DEFINITION_HPP
#define TABLE_DEFINITION_HPP

struct mytriplet {
  int i;
  int j;
  double v;
};

struct mytable {
  int Ti;
  int Tj;
  std::vector<mytriplet> tri;
};

#endif /* TABLE_DEFINITION_HPP */
