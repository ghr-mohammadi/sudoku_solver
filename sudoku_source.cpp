#include "sudoku_header.h"

int intRand(const int &min, const int &max) {
    static thread_local std::mt19937 *generator = nullptr;
    if (!generator)
        generator = new std::mt19937(clock() + (int) std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(*generator);
}

void printFunc(int *t)
{
	for (int j = 0; j < N; j++)
	{
		if (j % n == 0)
		{
			for (int k = 0; k < n * N - 2; k++)
				std::cout << "-";
			std::cout << std::endl;
		}
		for (int i = 0; i < N; i++)
		{
			if (i % n == 0)
				std::cout << "| ";
			if (t[i + j * N] > 0)
				std::cout << t[i + j * N] << " ";
			else
				std::cout << "* ";
		}
		std::cout << "|" << std::endl;
	}
	for (int k = 0; k < n * N - 2; k++)
		std::cout << "-";
	std::cout << std::endl;
}

int compare_(const void *a, const void *b) {
    return (*(int *) b - *(int *) a);
}

int init(int *t, int *t_) {
    int test[N], k;
    int i_test[N], j_test[N], pre_[N + 1];
    int i, j;
    int i_ = 0, j_ = 0;

    *(pre_ + N) = 0;

    for (int w = 0; w < N; w++)
        *(pre_ + w) = w;

    for (int _j = 0; _j < N; _j++) {
        int r = intRand(0, (N - 1) - _j);
        j_test[_j] = pre_[r];
        for (int kam = r; kam < N; kam++)
            pre_[kam] = pre_[kam + 1];
    }

    for (int j_it = 0; j_it < N; j_it++) {
        for (int w = 0; w < N; w++)
            *(pre_ + w) = w;

        for (int _i = 0; _i < N; _i++) {
            int r = intRand(0, (N - 1) - _i);
            i_test[_i] = pre_[r];
            for (int kam = r; kam < N; kam++)
                pre_[kam] = pre_[kam + 1];
        }

        for (int i_it = 0; i_it < N; i_it++) {
            j = *(j_test + j_it);
            i = *(i_test + i_it);
            if (!*(t_ + (i + j * N))) {
                for (int p = 0; p < N; p++)
                    *(test + p) = p + 1;

                for (int q = 0; q < N; q++)
                    if (*(t + (i + q * N)))
                        for (int p = 0; p < N; p++)
                            if (*(t + (i + q * N)) == *(test + p)) {
                                *(test + p) = 0;
                                break;
                            }

                for (int q = 0; q < N; q++)
                    if (*(t_ + (q + j * N)))
                        for (int p = 0; p < N; p++)
                            if (*(t_ + (q + j * N)) == *(test + p)) {
                                *(test + p) = 0;
                                break;
                            }

                j_ = j / n;
                i_ = i / n;
                for (int q = (j_ * n); q < ((j_ + 1) * n); q++)
                    for (int p = (i_ * n); p < ((i_ + 1) * n); p++)
                        if (*(t_ + (p + q * N)))
                            for (int r = 0; r < N; r++)
                                if (*(t_ + (p + q * N)) == *(test + r)) {
                                    *(test + r) = 0;
                                    break;
                                }

                qsort(test, N, sizeof(int), compare_);

                for (k = 0; *(test + k); k++);

                if (k == 0) {
                    j_ = j / n;
                    for (int q = (j_ * n); q < N; q++)
                        for (int p = 0; p < N; p++)
                            *(t_ + (p + q * N)) = *(t + (p + q * N));
                    return 0;
                }

                *(t_ + (i + j * N)) = *(test + intRand(0, k - 1));
            }
        }
    }
    return 1;
}

void initFunc(int *t, int *t_) {
    for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++)
            *(t_ + (i + j * N)) = *(t + (i + j * N));

    while (!init(t, t_));
}

bool numTest(char ch)
{
	if (ch > 47 && ch < 58)
		return true;
	else
		return false;
}

bool lineTest(char * lchar, int * lint)
{
	int j = 0;
	for (int i = 0; i < 2 * N; i++)
	{
		if (numTest(lchar[i]))
		{
			lint[j] = lchar[i] - '0';
			j++;
		}
		else if (lchar[i] != 32 && lchar[i] != '\0')
		{
			for (int k = 0; k < N; k++)
				lint[k] = 0;
			std::cout << "Please Enter Correct Input..." << std::endl;
			return false;
		}
	}
	return true;
}

void getFunc(int *t)
{
	char * charLine = nullptr;
	int * intLine = nullptr;

    std::cout << std::endl;
    std::cout << "Please Enter Sudoku Number Line by Line." << std::endl;
    std::cout << "Please Separate Numbers by Press Space." << std::endl;
    std::cout << "Please Separate Lines by Press Enter." << std::endl;
    std::cout << "Please Enter Zero in Blank Space." << std::endl;
    std::cout << std::endl;

	for (int j = 0; j < N; j++)
	{
		charLine = new char[5 * N];
		intLine = new int[N];
		std::cin.getline(charLine, 5 * N);
		if (!lineTest(charLine, intLine))
			j--;
		else
			for (int i = 0; i < N; i++)
				t[i + j * N] = intLine[i];
		delete[] charLine;
		charLine = nullptr;
		delete[] intLine;
		intLine = nullptr;
	}
}

int fitness(int *t) {
    int num = 0, mismatch = 0;

    for (int i = 0; i < N; i++)
        for (int p = 1; p <= N; p++) {
            num = 0;
            for (int j = 0; j < N && num < 2; j++)
                if (*(t + (i + j * N)) == p)
                    num++;
            if (num != 1)
                mismatch++;
        }

    return (N2 - mismatch);
}

void preSolve(int *t) {
    Numbers nums[N2];

    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            nums[i + j * N]._num = *(t + (i + j * N));
            nums[i + j * N]._Create();
        }
    }

    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            if (nums[i + j * N]._num == 0) {
                for (int k = 1; k <= N; k++) {
                    bool flag = true;
                    if (flag) {
                        for (int p = 0; p < N; p++) {
                            if (nums[p + j * N]._num == k) {
                                nums[i + j * N]._num_array[k - 1] = 0;
                                flag = false;
                                break;
                            }
                        }
                    }
                    if (flag) {
                        for (int q = 0; q < N; q++) {
                            if (nums[i + q * N]._num == k) {
                                nums[i + j * N]._num_array[k - 1] = 0;
                                flag = false;
                                break;
                            }
                        }
                    }
                    if (flag) {
                        for (int q = (j / n) * n; q < ((j / n) * n + n); q++) {
                            for (int p = (i / n) * n; p < ((i / n) * n + n); p++) {
                                if (nums[p + q * N]._num == k) {
                                    nums[i + j * N]._num_array[k - 1] = 0;
                                    flag = false;
                                    break;
                                }
                            }
                            if (!flag) {
                                break;
                            }
                        }
                    }

                }
            } else if (nums[i + j * N]._num_array != nullptr) {
                delete nums[i + j * N]._num_array;
                nums[i + j * N]._num_array = nullptr;
            }
        }
    }

    for (int m = 1; m <= N; m++) {
        for (int j = 0; j < N; j++) {
            int h = 0;
            for (int i = 0; i < N; i++) {
                if (nums[i + j * N]._num == 0) {
                    if (nums[i + j * N]._num_array[m - 1] == m) {
                        h++;
                    }
                }
            }
            if (h == 1) {
                for (int i = 0; i < N; i++) {
                    if (nums[i + j * N]._num == 0) {
                        if (nums[i + j * N]._num_array[m - 1] == m) {
                            delete nums[i + j * N]._num_array;
                            nums[i + j * N]._num_array = nullptr;
                            nums[i + j * N]._num = m;
                            *(t + (i + j * N)) = m;
                            if (fitness(t) == N2) {
                                return;
                            }
                            preSolve(t);
                            return;
                        }
                    }
                }
            }
        }
        for (int i = 0; i < N; i++) {
            int v = 0;
            for (int j = 0; j < N; j++) {
                if (nums[i + j * N]._num == 0) {
                    if (nums[i + j * N]._num_array[m - 1] == m) {
                        v++;
                    }
                }
            }
            if (v == 1) {
                for (int j = 0; j < N; j++) {
                    if (nums[i + j * N]._num == 0) {
                        if (nums[i + j * N]._num_array[m - 1] == m) {
                            delete nums[i + j * N]._num_array;
                            nums[i + j * N]._num_array = nullptr;
                            nums[i + j * N]._num = m;
                            *(t + (i + j * N)) = m;
                            if (fitness(t) == N2) {
                                return;
                            }
                            preSolve(t);
                            return;
                        }
                    }
                }
            }
        }
        for (int q = 0; q < n; q++) {
            for (int p = 0; p < n; p++) {
                int sq = 0, i = 0, j = 0;
                for (int j_ = 0; j_ < n; j_++) {
                    for (int i_ = 0; i_ < n; i_++) {
                        j = j_ + q * n;
                        i = i_ + p * n;
                        if (nums[i + j * N]._num == 0) {
                            if (nums[i + j * N]._num_array[m - 1] == m) {
                                sq++;
                            }
                        }
                    }
                }
                if (sq == 1) {
                    for (int j_ = 0; j_ < n; j_++) {
                        for (int i_ = 0; i_ < n; i_++) {
                            j = j_ + q * n;
                            i = i_ + p * n;
                            if (nums[i + j * N]._num == 0) {
                                if (nums[i + j * N]._num_array[m - 1] == m) {
                                    delete nums[i + j * N]._num_array;
                                    nums[i + j * N]._num_array = nullptr;
                                    nums[i + j * N]._num = m;
                                    *(t + (i + j * N)) = m;
                                    if (fitness(t) == N2) {
                                        return;
                                    }
                                    preSolve(t);
                                    return;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

int *solveFunc(int *t) {
    int t_[NP * N2], tm = 0;
    int swap_v = 0, changeFlag = 0, counter = 100, initFlag = 0;
    int a[NP], b[NP], c[NP];
    int good[(NP / 2) + 1], bad[(NP / 2) + 1];
    int v_[NP * N2], v_sort[NP * N2], temp_v[n * N];
    int temp_test[phi * N2], phi_counter = 0, temp_[N2];

    preSolve(t);

    if (fitness(t) == N2) {
        return (t);
    }

    while (true) {
        if (!(initFlag % tir)) {
            #pragma omp parallel for shared(t, t_) private(tm)
                for (tm = 0; tm < NP; tm++) {
                    initFunc(t, t_ + (tm * N2));
                }
            #pragma omp parallel for collapse(2) shared(temp_)
                for (int j = 0; j < N; j++) {
                    for (int i = 0; i < N; i++) {
                        *(temp_ + (i + (j * N))) = 0;
                    }
                }
        } else {
            #pragma omp parallel for shared(temp_, t_) private(tm)
                for (tm = 0; tm < NP; tm++) {
                    initFunc(temp_, t_ + (tm * N2));
                }
        }

        while (phi_counter < phi) {
            changeFlag = 0;

            for (int i = 0; i < NP; i++) {
                *(a + i) = fitness(t_ + (i * N2));
            }

            int crossOverTester = 0;
            for (int i = 1; i < NP; i++) {
                if ((*(a + i - 1)) != (*(a + i))) {
                    crossOverTester = 1;
                    break;
                }
            }

            if (crossOverTester) {
                for (int i = 0; i < NP; i++) {
                    *(b + i) = *(a + i);
                }

                qsort(b, NP, sizeof(int), compare_);

                *c = 1;
                for (int i = 1; i < NP; i++) {
                    if (!((*(b + i)) - *(b + i - 1))) {
                        *(c + i) = 0;
                    } else {
                        *(c + i) = 1;
                    }
                }

                for (int q = 0; q < NP; q++) {
                    int checker = 0;
                    for (int p = 0; p < NP; p++) {
                        if ((*(b + q)) == (*(a + p)) && *(c + q)) {
                            for (int j = 0; j < N; j++) {
                                for (int i = 0; i < N; i++) {
                                    *(v_sort + (i + (j * N) + ((q + checker) * N2))) = *(t_ + (i + (j * N) + (p * N2)));
                                }
                            }
                            checker++;
                        }
                    }
                }

                for (int k = 0; k < (NP - (int) (sNP * NP)); k++)
                    for (int j = 0; j < N; j++)
                        for (int i = 0; i < N; i++)
                            *(v_ + (i + (j * N) + (k * N2))) = *(v_sort + (i + (j * N) + (0 * N2)));

                for (int k = (NP - (int) (sNP * NP)); k < NP; k++)
                    for (int j = 0; j < N; j++)
                        for (int i = 0; i < N; i++)
                            *(v_ + (i + (j * N) + (k * N2))) = *(v_sort + (i + (j * N) + ((k - (NP - (int) (sNP * NP))) * N2)));

                for (int i = 0; i < NP / 2; i++) {
                    *(good + i) = i + 1;
                    *(bad + i) = i + 1 + NP / 2;
                }

                *(good + (NP / 2)) = *(bad + (NP / 2)) = 0;

                for (int q = 0; q < NP / 2; q++) {
                    int index_g = intRand(0, ((NP / 2) - 1) - q);
                    int index_b = intRand(0, ((NP / 2) - 1) - q);

                    int kg = (*(good + index_g) - 1);
                    int kb = (*(bad + index_b) - 1);

                    for (int index = index_g; index < NP / 2; index++)
                        *(good + index) = *(good + index + 1);
                    for (int index = index_b; index < NP / 2; index++)
                        *(bad + index) = *(bad + index + 1);

                    for (int p = 0; p < n; p++)
                        if (intRand(0, 1))
                            for (int j = p * n; j < (p + 1) * n; j++)
                                for (int i = 0; i < N; i++) {
                                    swap_v = *(v_ + (i + (j * N) + (kg * N2)));
                                    *(v_ + (i + (j * N) + (kg * N2))) = *(v_ + (i + (j * N) + (kb * N2)));
                                    *(v_ + (i + (j * N) + (kb * N2))) = swap_v;
                                }
                }

                for (int k = 0; k < NP; k++)
                    if (fitness(t_ + (k * N2)) < fitness(v_ + (k * N2))) {
                        for (int j = 0; j < N; j++)
                            for (int i = 0; i < N; i++)
                                *(t_ + (i + (j * N) + (k * N2))) = *(v_ + (i + (j * N) + (k * N2)));
                        changeFlag = 1;
                    } else
                        for (int j = 0; j < N; j++)
                            for (int i = 0; i < N; i++)
                                *(v_ + (i + (j * N) + (k * N2))) = *(t_ + (i + (j * N) + (k * N2)));

                for (int i = 0; i < NP; i++)
                    a[i] = fitness(v_ + (i * N2));
            }
            
            int blablaTester = 0;
            for (int i = 1; i < NP; i++)
                if ((*(a + i - 1)) != (*(a + i))) {
                    blablaTester = 1;
                    break;
                }

            int tedad = 1;
            if (!blablaTester) {
                if ((1.0 * (*a) / N2) > 0.95)
                    tedad = n * N2;
                else if ((1.0 * (*a) / N2) > 0.90)
                    tedad = N2;
                else if ((1.0 * (*a) / N2) > 0.85)
                    tedad = n * N;
                else if ((1.0 * (*a) / N2) > 0.80)
                    tedad = N;
            }

            if (blablaTester) {
                for (int u = 0; u < NP; u++)
                    for (int q = 0; q < n; q++) {
                        for (int j = q * n; j < (q + 1) * n; j++)
                            for (int i = 0; i < N; i++)
                                *(temp_v + (i + ((j - (q * n)) * N))) = *(v_ + (i + (j * N) + (u * N2)));

                        int test_i = 0, test_j = 0;

                        while (!(test_i - test_j)) {
                            test_i = intRand(1, N);
                            test_j = intRand(1, N);
                        }

                        for (int j = 0; j < n; j++)
                            for (int i = 0; i < N; i++)
                                if (*(temp_v + (i + (j * N))) == test_i)
                                    *(temp_v + (i + (j * N))) = test_j;
                                else if (*(temp_v + (i + (j * N))) == test_j)
                                    *(temp_v + (i + (j * N))) = test_i;

                        int tester = 0;
                        for (int k = 0; k < n; k++)
                            for (int num = 1; num <= N; num++) {
                                int checker = 0;
                                for (int j = 0; j < n; j++)
                                    for (int i = k * n; i < (k + 1) * n; i++)
                                        if (num == *(temp_v + (i + (j * N))))
                                            checker++;
                                if (checker != 1)
                                    tester++;
                            }

                        if (!tester) {
                            int flag = 0;

                            for (int j = q * n; j < (q + 1) * n; j++) {
                                for (int i = 0; i < N; i++)
                                    if ((*(t + (i + (j * N)))) * (*(temp_v + (i + ((j - (q * n)) * N))) - *(t + (i + (j * N))))) {
                                        flag = 1;
                                        break;
                                    }
                                if (flag)
                                    break;
                            }

                            if (!flag)
                                for (int j = q * n; j < (q + 1) * n; j++)
                                    for (int i = 0; i < N; i++)
                                        *(v_ + (i + (j * N) + (u * N2))) = *(temp_v + (i + ((j - (q * n)) * N)));
                        }
                    }

                for (int i = 0; i < NP; i++)
                    a[i] = fitness(v_ + (i * N2));

                for (int k = 0; k < NP; k++)
                    if (fitness(t_ + (k * N2)) < fitness(v_ + (k * N2))) {
                        for (int j = 0; j < N; j++)
                            for (int i = 0; i < N; i++)
                                *(t_ + (i + (j * N) + (k * N2))) = *(v_ + (i + (j * N) + (k * N2)));
                        changeFlag = 1;
                    } else
                        for (int j = 0; j < N; j++)
                            for (int i = 0; i < N; i++)
                                *(v_ + (i + (j * N) + (k * N2))) = *(t_ + (i + (j * N) + (k * N2)));
            } else {
                for (int u = 0; u < (tedad * N); u++)
                    for (int q = 0; q < n; q++) {
                        for (int j = q * n; j < (q + 1) * n; j++)
                            for (int i = 0; i < N; i++)
                                *(temp_v + (i + ((j - (q * n)) * N))) = *(v_ + (i + (j * N) + (0 * N2)));

                        int test_i = 0, test_j = 0;

                        while (!(test_i - test_j)) {
                            test_i = intRand(1, N);
                            test_j = intRand(1, N);
                        }

                        for (int j = 0; j < n; j++)
                            for (int i = 0; i < N; i++)
                                if (*(temp_v + (i + (j * N))) == test_i)
                                    *(temp_v + (i + (j * N))) = test_j;
                                else if (*(temp_v + (i + (j * N))) == test_j)
                                    *(temp_v + (i + (j * N))) = test_i;

                        int tester = 0;
                        for (int k = 0; k < n; k++)
                            for (int num = 1; num <= N; num++) {
                                int checker = 0;
                                for (int j = 0; j < n; j++)
                                    for (int i = k * n; i < (k + 1) * n; i++)
                                        if (num == *(temp_v + (i + (j * N))))
                                            checker++;
                                if (checker != 1)
                                    tester++;
                            }

                        if (!tester) {
                            int flag = 0;

                            for (int j = q * n; j < (q + 1) * n; j++) {
                                for (int i = 0; i < N; i++)
                                    if ((*(t + (i + (j * N)))) * (*(temp_v + (i + ((j - (q * n)) * N))) - *(t + (i + (j * N))))) {
                                        flag = 1;
                                        break;
                                    }
                                if (flag)
                                    break;
                            }

                            if (!flag)
                                for (int j = q * n; j < (q + 1) * n; j++)
                                    for (int i = 0; i < N; i++) {
                                        *(v_ + (i + (j * N) + (0 * N2))) = *(temp_v + (i + ((j - (q * n)) * N)));
                                    }
                        }
                    }

                if (fitness(t_ + (0 * N2)) < fitness(v_ + (0 * N2))) {
                    for (int j = 0; j < N; j++)
                        for (int i = 0; i < N; i++)
                            *(t_ + (i + (j * N) + (0 * N2))) = *(v_ + (i + (j * N) + (0 * N2)));
                    changeFlag = 1;
                    break;
                } else
                    for (int j = 0; j < N; j++)
                        for (int i = 0; i < N; i++)
                            *(v_ + (i + (j * N) + (0 * N2))) = *(t_ + (i + (j * N) + (0 * N2)));
            }

            counter--;

            for (int i = 0; i < NP; i++)
                if (fitness(t_ + (i * N2)) == N2) {
                    return (t_ + (i * N2));
                }

            if (counter == 0) {

                blablaTester = 0;
                for (int i = 1; i < NP; i++)
                    if ((*(a + i - 1)) != (*(a + i))) {
                        blablaTester = 1;
                        break;
                    }

                if (!blablaTester && (1.0 * (*a) / N2) > deghat) {
                    for (int j = 0; j < N; j++)
                        for (int i = 0; i < N; i++)
                            *(temp_test + (i + (j * N) + (phi_counter * N2))) = *(v_ + (i + (j * N) + (0 * N2)));
                    phi_counter++;
                }

                if (!(initFlag % tir)) {
                    #pragma omp parallel for shared(t, t_) private(tm)
                        for (tm = 0; tm < NP; tm++) {
                            initFunc(t, t_ + (tm * N2));
                        }
                } else {
                    #pragma omp parallel for shared(temp_, t_) private(tm)
                        for (tm = 0; tm < NP; tm++) {
                            initFunc(temp_, t_ + (tm * N2));
                        }
                }

                counter = 100;
            }
        }

        for (int j = 0; j < N; j++)
            for (int i = 0; i < N; i++) {
                int test_flag = 0;
                for (int k = 1; k < phi; k++)
                    if (*(temp_test + (i + (j * N) + ((k - 1) * N2))) !=
                        *(temp_test + (i + (j * N) + (k * N2)))) {
                        test_flag = 1;
                        break;
                    }
                if (!test_flag)
                    *(temp_ + (i + (j * N))) = *(temp_test + (i + (j * N) + (0 * N2)));
            }
        phi_counter = 0;
        initFlag++;
    }

}
