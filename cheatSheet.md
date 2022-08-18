find . -name 'error.dat'       -exec echo {} \; -exec tail -1 {} \; >> errors.txt;
find . -name 'convergence.dat' -exec echo {} \; -exec tail -1 {} \; >> convergence.txt;

watch -t "find . -name 'error.dat'       -exec "echo -n {}" \; -exec "tail -1 {}" \; | cut -c1-130 | sort"
watch -t "find . -name 'convergence.dat' -exec "echo -n {}" \; -exec "tail -1 {}" \; | cut -c1-100 | sort"
