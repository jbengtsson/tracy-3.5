#ifndef SIGMA_TRACK_H
#define SIGMA_TRACK_H

class bunch_train_type {
 private:

 public:
  int                         n_bunch;
  std::vector< ss_vect<tps> > bunch_train;
};

#endif
