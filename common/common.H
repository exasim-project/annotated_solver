#define START(NAME) auto start_##NAME = std::chrono::steady_clock::now()

#define STOP(NAME)                                                       \
	auto end_##NAME = std::chrono::steady_clock::now();              \
	Info << "[INFO] " #NAME ": " <<                                  \
	std::chrono::duration_cast<std::chrono::microseconds>            \
	(end_##NAME-start_##NAME).count()/1000.0 << " [ms]" << endl 
