#include <iostream>
#include <vector>
using namespace std;

template<typename T>

vector<float> linspace(T start_in, T end_in, int num_in)
{

  vector<float> linspaced;

  float start_ = static_cast<float>(start_in);
  float end_ = static_cast<float>(end_in);
  float num_ = static_cast<float>(num_in);

  if (num_ == 0) { return linspaced; }
  if (num_ == 1)
    {
      linspaced.push_back(start_);
      return linspaced;
    }

  float delta = (end_ - start_) / (num_ - 1);

  for(int i=0; i < num_-1; ++i)
    {
      linspaced.push_back(start_ + delta * i);
    }
  linspaced.push_back(end_);

  return linspaced;
}

template<typename T>
vector<float> linspace_interval(T start_in, T end_in, float del_in,int num_in)
{

  vector<float> linspaced;

  float start_ = static_cast<float>(start_in);
  float end_ = static_cast<float>(end_in);
  float del_ = static_cast<float>(del_in);
  int num_ = static_cast<int> (num_in);

  if (del_ == 0) { return linspaced; }
  int i = 0;
  for(i=0;i<num_;i++)
  {
      linspaced.push_back(start_ + del_ * i);
  }
  /*
  while(start_ + del_*i <= end_)
  {
      linspaced.push_back(start_ + del_ * i);
      i+=1 ;
  }*/
  return linspaced;
}

