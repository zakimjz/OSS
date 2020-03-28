//! \file uniform_freq.cpp - main() function
#include <iostream>
#include "database.h"
#include "uniform_freq_random_walk.h"

using namespace std;

char* datafile;
int minsup;
int uniq_pat_count;
int max_iter;

typedef ExPattern<int, int> PAT;
template<> vector<int> Database<PAT>::_no_data = Database<PAT>::set_static_data();    
template<> PatternFactory<PAT>* PatternFactory<PAT>::_instance = PatternFactory<PAT>::set_static_data();
//template<> PatternFactory<PAT>::set_static_data();

void print_usage(char *prog) {
  cerr<<"Usage: "<<prog<<" -d data-file -c count -s minsup"<<endl;
  exit(0);
}
/*! Parsing arguments */
void parse_args(int argc, char* argv[]) {
  if(argc<7) {
    print_usage(argv[0]);
  }

  for (int i=1; i < argc; ++i){
    if (strcmp(argv[i], "-d") == 0){
      datafile=argv[++i];
    }
    else if (strcmp(argv[i], "-s") == 0){
      minsup=atoi(argv[++i]);
      cout << "Minimum Support:" << minsup << " " << endl;
    }
    else if(strcmp(argv[i],"-c") == 0){
      max_iter=atoi(argv[++i]);
    }
    else{
      print_usage(argv[0]);
    }
  }
}//end parse_args()

/*! main function */
int main(int argc, char *argv[]) {

	bool zero_neighbors;
  parse_args(argc, argv);
  Database<PAT>* database;
  /* creating database and loading data */
  try {
    database = new Database<PAT>(datafile);
    database->set_minsup(minsup);
  }
  catch (exception& e) {
    
    cout << e.what() << endl;
  }
  database->remove_infrequent_edges();
	//database->print_database();
  /* creating random_walk_manager and starting walk */
  Uniform_Freq_Random_Walk<PAT> rwm(database, max_iter);
  int cur_iter=0;
	
	//will call initialize followed by walk again and again until it gets a single edge pattern with some neighbors.
	do {
			lattice_node<PAT>* start = rwm.initialize();      
			zero_neighbors = rwm.walk(start,cur_iter);
    } while (zero_neighbors == 0);
  delete database;
}
