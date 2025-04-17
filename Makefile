# Run instructions : "make 1" or "make 8" (any number from 1 to 8)

# # Suppress all output by default
.SILENT:

# Default target
all:
	@echo "Usage: make <1-8>"
	@echo "Example: make 1"

example:
	@if [ -f "src/example.cpp" ]; then \
		clear; \
		if $(MAKE) --no-print-directory --silent -C build example; then \
			./build/example; \
		else \
			echo "\n-------------------------------------------- \n\tBuild failed. Exiting...\n--------------------------------------------\n"; \
			exit 1; \
		fi; \
	else \
		echo "Error: example.cpp not found!"; \
		exit 1; \
	fi


# Pattern rule to handle targets 1-6
%:
	@if [ $(MAKECMDGOALS) -ge 1 ] && [ $(MAKECMDGOALS) -le 9 ]; then \
		clear; \
		if $(MAKE) --no-print-directory --silent -C build e$(MAKECMDGOALS); then \
			./build/e$(MAKECMDGOALS); \
		else \
			echo "\n-------------------------------------------- \n\tBuild failed. Exiting...\n--------------------------------------------\n"; \
			exit 1; \
		fi;\
	else \
		echo "Invalid target. Please provide a number between 1 and 8."; \
	fi