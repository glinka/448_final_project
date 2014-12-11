#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include "dmaps.h"
#include "gaussian_kernel.h"
#include "util_fns.h"
#include "embeddings.h"
#include "votingModel.h"
// #include "calcGraphProps.h"

int main(int argc, char** argv) {

  typedef std::vector< std::vector<int> > adj_mat;

  std::cout << "<---------------------------------------->" << std::endl;

  const int graph_size = 200;
  // run each model for total of total_nsteps steps, sampling every interval_nsteps steps for a total of nsamples graphs
  const int total_nsteps = atoi(argv[1]); // graph_size*graph_size;
  const int nsamples = atoi(argv[2]); //1000
  const int interval_nsteps = total_nsteps/nsamples;

  const double avg_deg = 4;
  const double a = 0.5;
  const int nopns = 2;
  const std::string rewireto = "random";
  double init_dist[2];
 
  // collect graphs from 5 different initial minority fractions (0.1 -> 0.5)
  const int ninits = 10;
  std::vector<adj_mat> vm_graphs(ninits*nsamples);
  std::vector< std::vector<int> > opn_embeddings(ninits*nsamples);
  int ngraphs = 0;

  // TESTING
  std::vector<double> minority_fracs(ninits*nsamples);
  std::vector<double> conflicts(ninits*nsamples);
  // END TESTING

  for(int i = 1; i < ninits + 1; i++) {
    init_dist[0] = 0.4 + i/(10.0*ninits);
    init_dist[1] = 1 - init_dist[0];
    votingModel votingmodel(graph_size, avg_deg, nopns, init_dist, a, rewireto);
    votingmodel.initGraph();

    bool converged = false;
    for(int j = 0; j < nsamples; j++) {
      if(!converged) {
  	vm_graphs[ngraphs++] = votingmodel.run_nsteps(interval_nsteps, converged);
  	opn_embeddings[ngraphs-1]  = votingmodel.get_opns();
	minority_fracs[ngraphs-1] = votingmodel.getMinorityFraction();
	conflicts[ngraphs-1] = votingmodel.getConflicts();
      }
    }
  }


  // collect samples of runs started from same initial conditions
  // const int ninits = 50;
  // init_dist[0] = 0.5;
  // init_dist[1] = 1 - init_dist[0];

  // for(int i = 0; i < ninits; i++) {
  //   votingModel votingmodel(graph_size, avg_deg, nopns, init_dist, a, rewireto);
  //   votingmodel.initGraph();

  //   bool converged = false;
  //   for(int j = 0; j < nsamples; j++) {
  //     if(!converged) {
  // 	vm_graphs[ngraphs++] = votingmodel.run_nsteps(interval_nsteps, converged);
  // 	opn_embeddings[ngraphs-1]  = votingmodel.get_opns();
  // 	minority_fracs[ngraphs-1] = votingmodel.getMinorityFraction();
  // 	conflicts[ngraphs-1] = votingmodel.getConflicts();
  //     }
  //   }
  // }

  //TESTING
  std::ofstream mfs_out("./single_runs/mf_evos.csv");
  std::ofstream cs_out("./single_runs/conflict_evos.csv");
  for(int i = 0; i < ngraphs; i++) {
    mfs_out << minority_fracs[i] << std::endl;
    cs_out << conflicts[i] << std::endl;
  }
  mfs_out.close();
  cs_out.close();
  // END TESTING

  
  std::cout << "--> " << ngraphs << " graph samples generated" << std::endl;

  // use both a spectral embedding and an "embedding" into a vector of opnions as input data

  // make evenly spaced number of spectral params (maybe should be logarithmically evenly spaced?)
  const int n_spec_params = 100;
  std::vector<double> spectral_params(n_spec_params);
  const double spec_param_max = 0.01;
  const double spec_param_min = 0.0001;
  const double dspec_param = (spec_param_max - spec_param_min)/(n_spec_params - 1);

  std::vector< std::vector<double> > graph_embeddings(ngraphs);
  for(int i = 0; i < n_spec_params; i++) {
    spectral_params[i] = spec_param_min + i*dspec_param;
  }
  for(int i = 0; i < ngraphs; i++) {
    graph_embeddings[i] = spectral_embedding(vm_graphs[i], spectral_params);
  }

  std::string label = "fullgraphs";
  std::cout << "--> Graphs embedded" << std::endl;
  std::cout << "--> Graph embeddings saved in: ./embedding_data" << std::endl;
  std::ofstream output_ges("./embedding_data/" + label + "_graph_embeddings.csv");
  utils::save_matrix(graph_embeddings, output_ges);
  output_ges.close();

  // DMAPS it
  const int k = std::atoi(argv[3]);
  std::vector<double> eigvals(k);
  std::vector< std::vector<double> > eigvects(k);
  std::vector< std::vector<double> > W(k);
  double median = utils::get_median(utils::get_squared_distances(graph_embeddings));
  std::cout << "--> Median squared distance: " << median << std::endl;
  Gaussian_Kernel<double> gk_d(median);
  int dmaps_success = dmaps::map(graph_embeddings, gk_d, eigvals, eigvects, W, k, 1e-12);
  if(dmaps_success != 1) {
    std::cout << "dmaps encountered an error" << std::endl;
  }

  std::ofstream output_eigvals("./embedding_data/dmaps_" + label + "_embedding_eigvals.csv");
  std::ofstream output_eigvects("./embedding_data/dmaps_" + label + "_embedding_eigvects.csv");

  output_eigvals << "ninits=" << ninits << std::endl;
  output_eigvects << "ninits=" << ninits << std::endl;
  utils::save_vector(eigvals, output_eigvals);
  utils::save_matrix(eigvects, output_eigvects);
  output_eigvals.close();
  output_eigvects.close();

  std::cout << "--> First DMAP computed" << std::endl;
  std::cout << "--> First DMAP output saved in: ./embedding_data" << std::endl;

  label = "minority_vects";
  // DMAPS it
  // remove unused spaces in vector
  opn_embeddings.erase(opn_embeddings.begin() + ngraphs, opn_embeddings.end());
  median = utils::get_median(utils::get_squared_distances(opn_embeddings));
  std::cout << "--> Median squared distance: " << median << std::endl;
  Gaussian_Kernel<int> gk_i(median);
  dmaps_success = dmaps::map(opn_embeddings, gk_i, eigvals, eigvects, W, k, 1e-12);
  if(dmaps_success != 1) {
    std::cout << "dmaps encountered an error" << std::endl;
  }

  output_eigvals.open("./embedding_data/dmaps_" + label + "_embedding_eigvals.csv");
  output_eigvects.open("./embedding_data/dmaps_" + label + "_embedding_eigvects.csv");

  output_eigvals << "ninits=" << ninits << std::endl;
  output_eigvects << "ninits=" << ninits << std::endl;
  utils::save_vector(eigvals, output_eigvals);
  utils::save_matrix(eigvects, output_eigvects);
  output_eigvals.close();
  output_eigvects.close();

  std::cout << "--> Second DMAP computed" << std::endl;
  std::cout << "--> Second DMAP output saved in: ./embedding_data" << std::endl;

  // DMAPS is hard, PCA it
  Eigen::MatrixXd X(ngraphs, n_spec_params);
  for(int i = 0; i < ngraphs; i++) {
    for(int j = 0; j < n_spec_params; j++) {
      X(i,j) = graph_embeddings[i][j];
    }
  }
  Eigen::MatrixXd V;
  Eigen::VectorXd S;
  const int pca_success = eigen_solver::pca(X, V, S, k);
  if(pca_success != 1) {
    std::cout << "pca encountered an error" << std::endl;
    return 0;
  }
  std::cout << "--> PCA completed" << std::endl;
  std::cout << "--> PCA output saved in: ./embedding_data" << std::endl;

  output_eigvals.open("./embedding_data/pca_" + label + "_embedding_eigvals.csv");
  output_eigvects.open("./embedding_data/pca_" + label + "_embedding_eigvects.csv");

  utils::save_vector(std::vector<double>(S.data(), S.data() + k), output_eigvals);
  for(int i = 0; i < k; i++) {
    eigvects[i] = std::vector<double>(V.data() + i*n_spec_params, V.data() + (i+1)*n_spec_params);
  }
  utils::save_matrix(eigvects, output_eigvects);
  output_eigvals.close();
  output_eigvects.close();

  std::cout << "<---------------------------------------->" << std::endl;
  
  return 0;
}
