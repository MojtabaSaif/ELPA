	


#if !defined(MUTUAL_INCLUDED)
#define MUTUAL_INCLUDED	
#include "histograms.h"
namespace NMI {

 class Mutal{
 public:
	int overlap_grouping(deque<deque<int> > ten, int unique) {		//hrepiguhpueh


		set <int> conta;
		int all = 0;


		for (int i = 0; i < ten.size(); i++) {

			for (int j = 0; j < ten[i].size(); j++)
				conta.insert(ten[i][j]);

			all += ten[i].size();

		}


		unique = conta.size();

		int overlap = all - unique;




		return overlap;

	}



	double mutual(deque<deque<int> > en, deque<deque<int> > ten) {


		// en e ten are two partitions of integer numbers
		int dim;


		{
			set <int> conta;
			//set <int> ten_;
			//set <int> en_;

			for (int i = 0; i < ten.size(); i++) {
				sort(ten[i].begin(), ten[i].end());
				for (int j = 0; j < ten[i].size(); j++) {
					conta.insert(ten[i][j]);
					//ten_.insert(ten[i][j]);
				}
			}

			for (int i = 0; i < en.size(); i++) {
				sort(en[i].begin(), en[i].end());
				for (int j = 0; j < en[i].size(); j++) {
					conta.insert(en[i][j]);
					//en_.insert(en[i][j]);

				}
			}




			dim = conta.size();

			/*
			for (set<int>::iterator its=conta.begin(); its!=conta.end(); its++) {

				if(ten_.find(*its)==ten_.end()) {

					deque <int> first;
					first.push_back(*its);
					ten.push_back(first);
				}

				if(en_.find(*its)==en_.end()) {

					deque <int> first;
					first.push_back(*its);
					en.push_back(first);
				}


			}
			*/


		}


		//cout<<"dim:\t"<<dim<<endl;


		deque<deque<double> > N;
		deque <double> first;
		first.assign(en.size(), 0);
		for (int i = 0; i < ten.size(); i++)
			N.push_back(first);

		deque <int> s(dim);
		for (int i = 0; i < ten.size(); i++)
			for (int j = 0; j < en.size(); j++)
				N[i][j] = set_intersection(ten[i].begin(), ten[i].end(), en[j].begin(), en[j].end(), s.begin()) - s.begin();

		//printv(N);


		/*
		cout<<"one:"<<endl;
		printm(ten);


		cout<<"two:"<<endl;
		printm(en);


		cout<<"confusion matrix"<<endl;
		printm(N, cout);
		*/


		//cout<<"confusion matrix"<<endl;
		//printm(N, cout);







		deque <double> NR;
		NR.assign(ten.size(), 0);
		deque <double> NC;
		NC.assign(en.size(), 0);
		double NTOT = dim;

		for (int i = 0; i < ten.size(); i++)
			for (int j = 0; j < en.size(); j++) {
				NR[i] += N[i][j];
				NC[j] += N[i][j];
			}



		double IN_ = 0;
		double ID1 = 0;
		double ID2 = 0;

		for (int i = 0; i < ten.size(); i++)
			for (int j = 0; j < en.size(); j++)
				if (N[i][j] != 0)
					IN_ += N[i][j] * log(N[i][j] * NTOT / (NR[i] * NC[j]));

		IN_ = -2. * IN_;


		for (int i = 0; i < ten.size(); i++)
			if (NR[i] != 0)
				ID1 += NR[i] * log(NR[i] / (NTOT));



		for (int j = 0; j < en.size(); j++)
			if (NC[j] != 0)
				ID2 += NC[j] * log(NC[j] / (NTOT));




		double I = IN_ / (ID1 + ID2);

		if ((ID1 + ID2) == 0)
			I = -2;


		return I;

	}


	double H(double a) {


		if (a <= 0)
			return 0;
		else
			return (-a * log(a));

	}


	double H(deque <double>& p) {

		double h = 0;
		for (deque<double>::iterator it = p.begin(); it != p.end(); it++) if (*it != 0)
			h += (*it) * log(*it);

		return (-h);


	}



	double H_x_given_y3(deque<deque<int> >& en, deque<deque<int> >& ten, int dim) {

		// you know y and you want to find x according to a certain index labelling.
		// so, for each x you look for the best y.


		deque<deque<int> > mems;

		deque<int> first;
		for (int i = 0; i < dim; i++)
			mems.push_back(first);

		for (int ii = 0; ii < ten.size(); ii++)
			for (int i = 0; i < ten[ii].size(); i++)
				mems[ten[ii][i]].push_back(ii);



		double H_x_y = 0;
		double H2 = 0;

		for (int k = 0; k < en.size(); k++) {



			deque<int>& c = en[k];


			deque <double> p;
			double I2 = double(c.size());
			double O2 = (dim - I2);
			p.push_back(I2 / dim);
			p.push_back(O2 / dim);
			double H2_ = H(p);
			p.clear();



			double diff = H2_;

			// I need to know all the group wuth share nodes with en[k]

			map<int, int> com_ol;		// it maps the index of the ten into the overlap with en[k]
			Histograms hist;
			for (int i = 0; i < c.size(); i++) {

				for (int j = 0; j < mems[c[i]].size(); j++)
					hist.int_histogram(mems[c[i]][j], com_ol);

			}



			for (map<int, int>::iterator itm = com_ol.begin(); itm != com_ol.end(); itm++) {



				double I1 = double(ten[itm->first].size());
				double O1 = (dim - I1);


				p.push_back(I1 / dim);
				p.push_back(O1 / dim);
				double H1_ = H(p);
				p.clear();




				double I1_I2 = itm->second;
				double I1_02 = ten[itm->first].size() - I1_I2;
				double O1_I2 = c.size() - I1_I2;
				double O1_02 = dim - I1_I2 - I1_02 - O1_I2;


				p.push_back(I1_I2 / dim);
				p.push_back(O1_02 / dim);

				double H12_positive = H(p);

				p.clear();
				p.push_back(I1_02 / dim);
				p.push_back(O1_I2 / dim);


				double H12_negative = H(p);

				double H12_ = H12_negative + H12_positive;

				p.clear();

				if (H12_negative > H12_positive) {

					H12_ = H1_ + H2_;

				}





				if ((H12_ - H1_) < diff) {
					diff = (H12_ - H1_);
				}



			}



			if (H2_ == 0)
				H_x_y += 1;
			else
				H_x_y += (diff / H2_);

		}







		return (H_x_y / (en.size()));

	}


	double mutual3(deque<deque<int> > en, deque<deque<int> > ten) {



		if (en.size() == 0 || ten.size() == 0)
			return 0;



		// en e ten are two partitions of integer numbers
		int dim;



		{
			map <int, int> all;		// node, index 
			//set <int> ten_;
			//set <int> en_;

			for (int i = 0; i < ten.size(); i++) {
				for (int j = 0; j < ten[i].size(); j++) {
					all.insert(make_pair(ten[i][j], all.size()));
					//ten_.insert(ten[i][j]);
				}
			}

			for (int i = 0; i < en.size(); i++) {
				for (int j = 0; j < en[i].size(); j++) {
					all.insert(make_pair(en[i][j], all.size()));

				}
			}




			dim = all.size();



			for (int i = 0; i < ten.size(); i++) {
				for (int j = 0; j < ten[i].size(); j++)
					ten[i][j] = all[ten[i][j]];

				sort(ten[i].begin(), ten[i].end());
			}

			for (int i = 0; i < en.size(); i++) {
				for (int j = 0; j < en[i].size(); j++)
					en[i][j] = all[en[i][j]];

				sort(en[i].begin(), en[i].end());
			}





		}




		return (0.5 * (2. - H_x_given_y3(ten, en, dim) - H_x_given_y3(en, ten, dim)));

	}

	};
}

#endif

