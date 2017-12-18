template <typename T>
class BaseHistogram {

	std::map<T, size_t> sorted_elems_;
	std::vector<T> part_sums_;
	bool invalid_part_sums_;

	void CalculatePartSums(){
		part_sums_.clear();
		T prev_elem = 0;
		for(auto it = sorted_elems_.begin(); it != sorted_elems_.end(); it++){
			T key = it->first;
			size_t count = it->second;
			for(size_t i = 0; i < count; i++){
				part_sums_.push_back(key + prev_elem);
				prev_elem += key;
			}
		}
		invalid_part_sums_ = false;
	}

	T operator[](size_t idx){
		VERIFY(idx < part_sums_.size());
		if(size() == 0)
			return 0;
		if(idx == 0)
			return part_sums_[0];
		return part_sums_[idx] - part_sums_[idx - 1];
	}

public:
	void Add(T new_elem, size_t count = 1){
		invalid_part_sums_ = true;
		if(sorted_elems_.find(new_elem) == sorted_elems_.end())
			sorted_elems_[new_elem] = count;
		else
			sorted_elems_[new_elem] += count;
	}

	T Quantile(double quantile){
		VERIFY(quantile > 0 && quantile <= 1);
		if(invalid_part_sums_)
			CalculatePartSums();
		if(part_sums_.size() == 0)
			return 0;
		T total_sum = part_sums_[part_sums_.size() - 1];
		for(size_t i = 0; i < part_sums_.size(); i++)
			if(double(part_sums_[i]) / double(total_sum) >= quantile)
				return (*this)[i];

		return T(0);
	}

	size_t size() { return part_sums_.size(); }

	T max() {
		CalculatePartSums();
		if(size() == 0)
			return 0;
		return (*this)[size() - 1];
	}

	void SaveToFile(std::string filename) const {
		std::ofstream out(filename.c_str());
 		for(auto it = sorted_elems_.begin(); it != sorted_elems_.end(); it++)
			out << it->first << ' ' << it->second << endl;
	}

	void LoadFrom(std::string filename) {
		ifstream src(filename.c_str());
		VERIFY(!src.fail());
		while(!src.eof()){
			std::string tmp;
			getline(src, tmp);
			if(tmp != ""){
				std::stringstream ss;
				ss << tmp;
				T elem;
				size_t count;
				ss >> elem;
				ss >> count;
				Add(elem, count);
			}
		}
	}

	void SimpleSavingToFile(std::string filename) const {
		std::ofstream out(filename.c_str());
 		for(auto it = sorted_elems_.begin(); it != sorted_elems_.end(); it++)
 			for(size_t i = 0; i < it->second; i++)
 				out << it->first << endl;
	}
};
