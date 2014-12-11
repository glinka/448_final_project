namespace utils {


  template <typename T>
  double average(const std::vector<T>& v) {
    double average = 0;
    int n = v.size();
    for(int i = 0; i < n; i++) {
      average += v[i];
    }
    return average/n;
  }

  template <typename T>
  double get_median(const std::vector<T>& v) {
    std::vector<double> vcopy = v;
    std::sort(vcopy.begin(), vcopy.end());
    int n = vcopy.size();
    if(n % 2 == 0) {
      return (vcopy[n/2] + vcopy[n/2 + 1])/2.0;
    }
    else {
      return vcopy[n/2 + 1];
    }
  }

  template <typename T>
  std::vector<double> get_squared_distances(const std::vector< std::vector<T> >& vectors) {
    const int nvects = vectors.size();
    int ncombos = (nvects*(nvects-1))/2;
    std::vector<double> squared_distances(ncombos);
    int counter = 0;
    for(int i = 0; i < nvects; i++) {
      for(int j = i+1; j < nvects; j++) {
	squared_distances[counter++] = std::pow(l2_norm(vectors[i], vectors[j]), 2);
      }
    }
    return squared_distances;
  }
   
  template <typename T>
  double l2_norm(const std::vector<T>& x1, const std::vector<T>& x2) {
    double norm = 0;
    for(int i = 0; i < x1.size(); i++) {
      norm += pow(x1[i] - x2[i], 2);
    }
    return sqrt(norm);
  }

  template <typename T>
  void save_matrix(const std::vector< std::vector<T> >& A, std::ofstream& output_file, const std::string header, const char delim) {
    if(!header.empty()) {
      output_file << header << std::endl;
    }
    for(typename std::vector< std::vector<T> >::const_iterator v = A.begin(); v != A.end(); v++) {
      save_vector(*v, output_file, "", delim);
    }
  }

  template <typename T>
  void save_vector(const std::vector<T>& v, std::ofstream& output_file, const std::string header, const char delim) {
    if(!header.empty()) {
      output_file << header << std::endl;
    }

    for(typename std::vector<T>::const_iterator val = v.begin(); val != v.end()-1; val++) {
      output_file << *val << delim;
    }
    // always ends with newline
    output_file << v.back() << std::endl;
  }  

}
