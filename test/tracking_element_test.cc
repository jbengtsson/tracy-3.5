// Element-level tracking regression tests.  Each group of cases lives in its
// own tracking_*.cc file; tracking_support.h holds what they share.
#include "tracking_support.h"

#include <cstdio>
#include <cstdlib>

int main()
{
  using namespace tracking_test;

  configure_tracking();
  run_element_tests();

  if (failure_count != 0) {
    std::fprintf(stderr, "%d tracking element regression failure(s)\n",
                 failure_count);
    return EXIT_FAILURE;
  }

  std::printf("tracking element regressions passed\n");
  return EXIT_SUCCESS;
}
