#include <string>
#include <map>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include <ax_map_analyser.hpp>
#include <ax_map_generator.hpp>

#define ANALISYS 1 // 1 = yes, 0 = nope
#define GENERATE_MAP 0

template<int N = 10>
class Hystory {
private:
int array_[N] = {0};
int idx_ = 0;

public:
int sum() {
    int res = 0;
    for (int i = 0; i < N; i++)
        res += array_[i];
    return res;
}

void push(int x) {
    array_[idx_] = x;
    idx_ = (idx_ + 1) % N;
}
};

using namespace ax;

float mapD04(int x, int y)
{
    int X = 20, Y = 20;
    float R_SAN = 6, R_AT = 8;
    float sigma = 0.6;
    float r = r_points(Point(X / 2, Y / 2), Point(x, y));
    float type;

    if (r <= R_SAN)
    { // inside
        do {
            type = Get_Normal(r / R_SAN, sigma);
        }  while (type < 0 || type > 1);
        return type;
    } else if (r <= R_AT) { // buffer zone
            float reverse_r = 1.f - (r - R_SAN) / (R_AT - R_SAN);
            if (fabsf(Get_Normal<float>(0, 0.7f)) > reverse_r)
                return mag::Cell_AT();
            do {
                type = Get_Normal<float>(1, sigma);
            } while (type < 0 || type > 1);
            return type;
        }
    return mag::Cell_AT();
}

int main(int argc, char *argv[]) {

#if GENERATE_MAP == 1
    // creating map

  // added seed
  int seed = 8;
  srand(seed);
  mag::Map tmap(20, 20, "map_200x200_karpaev.map");
  tmap.Create(mapD04);
  //tmap.Filter(ax::FILTER_sharp);
  tmap.Write();

  std::fstream forSeed("seed.txt", std::ios::out);
  forSeed << "Seed = " << seed << std::endl;
#endif


#if ANALISYS == 1
    // pmaking analisys
    /*man::Interpolant2 I2(0, 1, 0);
	printf("%f\n", I2(-1));
	printf("%f\n", I2.Reverse(0.75));
  
	man::Interpolant1 I1(0, 1);
	printf("%f\n", I1(-1));
	printf("%f\n", I1.Reverse(0.5));
  
	printf("\n");
	man::Interpolant3 I3(0, 1, -1, 0);
	printf("%f\n", I3(-2.f));
	printf("%f\n", I3(-1.5f));
	printf("%f\n", I3(-1.f));
	printf("%f\n", I3.Reverse(-0.5f));
	  
	*/

    if (argc < 4)
        return 1;
    int X, Y;
    X = atoi(argv[1]);
    Y = atoi(argv[2]);

    man::Map<10> map(argv[3], X, Y);
    std::string oname(argv[3]);

    std::string fileFormat = ".txt", timePeriodString = "_less_1sec";
    oname += timePeriodString + fileFormat;
    FILE *fout = fopen(oname.c_str(), "w");
    int edge = 3;

    Hystory<20> hystory;
    auto STATE = man::PROCESS_STATE::RELAXING;
    std::list<std::pair<int, int>> stacked_centres; // (x, y)
    std::map<std::pair<int, int>, int> all_centres; // (x, y, times)

    // -- my code
    std::array<std::array<int, 200>, 200> CurrentLCpresenceMatrix, PreviousLCpresenceMatrix;
    std::array<std::array<float, 200>, 200> LCpresenceWholeTimeMatrix;

    for (auto &rows: LCpresenceWholeTimeMatrix)
        for (auto &elem: rows)
            elem = 0.;

    for (auto &rows: CurrentLCpresenceMatrix)
        for (auto &elem: rows)
        {
            elem = 0;
            std::cout << "elem = " << elem << std::endl;
        }

    // /exit(1);
    std::vector<std::tuple<float, int, int>> all_centres_repeat; // (t, x, y)



    Point corner(7, 7);

    //std::list<std::pair<int, int>> V_points = {{60, 100}, {115, 100}, {100, 85}, {100, 140}};

    bool success = true;     /// TOCHANGE if analysing before/after ACH period!
    Moving_average T(0/*0.5f*/);
    std::vector<float> listOfCycleLengths;
    int numCycles = 0;

    fprintf(fout, "время, с||кол-во ВЦ||период колебаний, с\n");
    map.ClearOutFrame();
    float *tempOut = new float[X * Y];
    double tStop = 1., tStart = 0;

    float last_T = tStart, first_T;
    
    bool isLC = false;



    // timeloop
    for (int i = 0; success; i++)
    {
        
        if (i > tStop/ map.GetTimestep())
            break;

        if (i < tStart/map.GetTimestep())
        {	
        	if (i % 100 == 0)
				std::cout << "Skipping frame #" << i << std::endl;
        	success = map.NextFrame();
        	continue;
        }

        decltype(stacked_centres) current_centres;

        /*if (i <  60./ map.GetTimestep())
		{
		  if (i % 1000 == 0)
				std::cout << "Skipping frame #" << i << std::endl;
			success = map.NextFrame();
			continue;
		}
		
		if (i > 24000)
		  break;*/


        /*/! V
		float V;
		for (auto &&p : V_points) {
		  float g_mod;
		  auto grad = map.getGradient(p.first, p.second, &g_mod);
		  V = map.getWaveSpeed(p.first, p.second);
		  if (V > 0)   printf("\n%d (%d %d): grad = (%f %f), V = %f\n",
				   map.getFrame(), p.first, p.second, grad.first, grad.second, V);
		}
		//*/
        map.ClearOutFrame(); // LC
        
        // going through all cells
        for (int y = edge; y < Y - edge; ++y) 
        {
            for (int x = edge; x < X - edge; ++x) 
            {
                //printf("x, y: (%d, %d)\n", x, y);
                //! Leading centre
                isLC = map.isLeadingCentreNew(x, y);
                if (isLC) 
                {
                    current_centres.emplace_back(x, y);
                    all_centres_repeat.emplace_back(std::make_tuple(map.getTime(), x, y));

                    int EoutOld = map.Eout(x, y);
                    map.Eout(x, y) += 1; // LC number at point (x, y)
                    CurrentLCpresenceMatrix[x][y] = 1;
                    //return -1;  // FOR DEBUG
                }

                //*/

                /*/! Grad
				float g;
				map.getGradient(x, y, &g);
				map.Eout(x, y) = g;
				 //*/

                //if (map.isSpark(x, y)) map.Eout(x, y) = map.getTime(); // chronotopographic map

                /*/! Reentry
				  map.Eout(x, y) = map.isSingularity(x, y);
				//*/

                //map.Eout(x, y) = map.getPhase(x, y); // Phase
            }
        }

        /*/! Grad
		for (int y = edge; y < Y - edge; ++y) {
		  for (int x = edge; x < X - edge; ++x) {
			float g;
			auto grad = map.getGradient(x, y, &g);
			if (x % 7 == 0 && y % 7 == 0) {
			  map.Write_Line_Out_Frame({x, y}, grad, {-80, 30}, (int)(g/2));
			}
	  
	  }
		}
		//*/
        //map.WriteOutFrame(); // Phase
        /*/! Reentry
		map.FilterOutFrame(man::FILTER_mean3x3);
		for (int y = edge; y < Y - edge; ++y) {
		  for (int x = edge; x < X - edge; ++x) {
			if (map.Eout(x, y) > 0.95f) {
			  fprintf(fout, "%d,%.6f,%d,%d,%.2f\n",
					  map.getFrame(), map.getTime(), x, y, map.Eout(x, y));
			}
		  }
		}
		fflush(fout);
		 */
        //map.WriteOutFrame();

        /*// Interpolants test
		for(float t = -1.0f; t < 0.0f; t += 0.1f) {
		  map.ClearOutFrame();
		  for (int y = 0; y < Y; ++y) {
			for (int x = 0; x < X; ++x) {
			  map.Eout(x, y) = map.interpols1_[map.xy(x, y)](t);
			}
		  }
		  map.WriteOutFrame();
		}*/


        hystory.push(static_cast<int>(current_centres.size()));
        
        if (current_centres.size() > 0)
        {
            /*if (STATE == man::PROCESS_STATE::RELAXING)
            {
                numCycles++;
                if (numCycles <= 1)
                {
                    last_T = map.getTime();  // starting time of calculation of LC denssity
                    startTime = last_T;

                    for (int x = 0; x < X; x++) // filling Prev. LC presence matrix for 1st time
                        for (int y = 0; y < Y; y++)
                            PreviousLCpresenceMatrix[x][y] = CurrentLCpresenceMatrix[x][y];

                    for (auto &rows: CurrentLCpresenceMatrix)
                        for (auto &elem: rows) 
                        {
                            elem = 0;
                        }

                }

                else 
                {
                    T.update(map.getTime() - last_T);   // calculate cycle length
                    listOfCycleLengths.push_back(T.get());

                    // filling for density map
                    for (int x = 0; x < X; x++)
                        for (int y = 0; y < Y; y++)
                            LCpresenceWholeTimeMatrix[x][y] += listOfCycleLengths.back() * PreviousLCpresenceMatrix[x][y];


                    for (int x = 0; x < X; x++) // changing LC presence layers
                        for (int y = 0; y < Y; y++)
                            PreviousLCpresenceMatrix[x][y] = CurrentLCpresenceMatrix[x][y];

                    
                    for (auto &rows: CurrentLCpresenceMatrix)
                        for (auto &elem: rows) {
                            elem = 0;
                        }

                }
            }*/

            

            /*for (auto &val: stacked_centres)
			{
				std::cout << "Num_centres = " << std::get<0>(val) << std::endl;
			}*/

            STATE = man::PROCESS_STATE::IGNITION;
        }

        stacked_centres.splice(stacked_centres.end(), current_centres);



        if (STATE == man::PROCESS_STATE::IGNITION && hystory.sum() == 0) 
        { // IGNITION -> RELAXING, drop stacked
            
            numCycles++;
            if (numCycles <= 1)
            {
                    last_T = map.getTime();  // starting time of calculation of LC denssity
                    first_T = last_T;

                    for (int x = 0; x < X; x++) // filling Prev. LC presence matrix for 1st time
                        for (int y = 0; y < Y; y++)
                            PreviousLCpresenceMatrix[x][y] = CurrentLCpresenceMatrix[x][y];

                    for (auto &rows: CurrentLCpresenceMatrix)
                        for (auto &elem: rows) 
                        {
                            elem = 0;
                        }

            } // comment if 1st ignition to be counted

            else
            {
                 T.update(map.getTime() - last_T);
                 last_T = map.getTime();
                // calculate cycle length
                 listOfCycleLengths.push_back(T.get());

                 // filling for density map
                 for (int x = 0; x < X; x++)
                     for (int y = 0; y < Y; y++)
                        LCpresenceWholeTimeMatrix[x][y] += listOfCycleLengths.back() * PreviousLCpresenceMatrix[x][y];


                for (int x = 0; x < X; x++) // changing LC presence layers
                    for (int y = 0; y < Y; y++)
                       PreviousLCpresenceMatrix[x][y] = CurrentLCpresenceMatrix[x][y];

                    
                 for (auto &rows: CurrentLCpresenceMatrix)
                    for (auto &elem: rows) {
                         elem = 0;
                }
            }







            fprintf(fout, "%.6f,%d,%.4f\n",
            /*map.getFrame(),*/
            map.getTime(),
            static_cast<int>(stacked_centres.size()),
            T.get());
            fflush(fout);

            for (auto &&c : stacked_centres) {
                all_centres[c]++;
            }

            stacked_centres.clear();
            //all_centres_repeat.splice(all_centres_repeat.end(), stacked_centres);
            STATE = man::PROCESS_STATE::RELAXING;


        




            //! Leading centres
            { // Frame with density
                //map.FilterOutFrame(ax::FILTER_mean3x3);
                map.NormalizeOutFrame(100);
                //map.WriteOutFrame();
                //map.Draw_Number_Out_Frame<100>(corner, map.getFrame());
                map.DumpOutFrame(tempOut);
                map.ClearOutFrame();
            } //*/
        }

        //! Leading centres
        if (i % 5000 == 0) // every 0.5 s
            map.WriteOutFrame(tempOut);
        //*/

        if (map.getFrame() % 10 == 0)
            printf("|");
        if (map.getFrame() % 1000 == 0)
            printf("\nframe: %d\n", map.getFrame());
        fflush(stdout);

        success = map.NextFrame();
    };

    //fclose(fout);

    // DEBUG
    //for (auto const& rows: LCpresenceWholeTimeMatrix)
    //for (auto const& elem: rows)
    //std::cout << "matrix = " << 100 * elem / 60. << '\n';




    // next section
    std::string oname2(argv[3]);
    oname2 += "all_centres_" + timePeriodString + fileFormat;
    FILE *fout2 = fopen(oname2.c_str(), "w");
    double threshold = (map.getTime() / T.get()) * 0.01; // 1% of all ignitions
    //fprintf(fout, "\n\n--------------------- %s ---------------------\n", "Stat: all centres with number (x, y, times)");
    //fprintf(fout, "Threshold = %d\n", static_cast<int>(threshold));
    fprintf(fout2, "x, клеток||y, клеток||плотность\n");
    // { plot fix
    fprintf(fout2, "%d,%d,%d\n", 0, 0, 0);
    fprintf(fout2, "%d,%d,%d\n", 0, Y, 0);
    fprintf(fout2, "%d,%d,%d\n", X, 0, 0);
    fprintf(fout2, "%d,%d,%d\n", X, Y, 0);
    // } fix
    std::list<std::pair<std::pair<int, int>, int>> list_centres;
    for (auto &&c : all_centres) {
        if (c.second < threshold)
            continue;
        list_centres.emplace_back(c.first, c.second);
    }
    list_centres.sort([](std::pair<std::pair<int, int>, int> &a, std::pair<std::pair<int, int>, int> &b) -> bool {
    return a.second > b.second;
    });

    map.ClearOutFrame();
    for (auto &&c : list_centres) {
        fprintf(fout2, "%d,%d,%d\n", c.first.first, c.first.second, c.second);
        map.Eout(c.first.first, c.first.second) = c.second;
    }
    //map.FilterOutFrame(ax::FILTER_mean3x3);
    map.NormalizeOutFrame(100);
    map.WriteOutFrame();

    //fclose(fout2);



    std::string fn(argv[3]);
    fn += "_LC_"+ timePeriodString + fileFormat;
    //FILE *fout = fopen(oname.c_str(), "w");
    //int edge = 10;
    // VTK output for Leading Centers
    //char fn[256];
    //sprintf(vt, "LeadingCenters.vtk");

    std::fstream fLC(fn, std::ios::out);
    double value;

    for (int x = 0; x < X; x++)
        for (int y = 0; y < Y; y++)
        {
            value = LCpresenceWholeTimeMatrix[x][y];
            fLC << x << "," << y << "," << value << "\n";  // "/100." -- back to fraction instead of %
        }
    fLC.close();

    /*
	for (int j = 0; j < atoi(argv[1]); j++) 
	  {
		  for (int i = 0; i < atoi(argv[2]); i++)
			  fLC << map.Eout(i, j) << " ";
	  }
	fLC << std::endl;
	fLC.close();*/

    /*
	fprintf(fout, "\nMAP\n");
	for (int x = 0; x < X; ++x)
	  fprintf(fout, "%d,", x + 1);
	fprintf(fout, "\n");
	for (int y = 0; y < Y; ++y) {
	  fprintf(fout, "%d,", y + 1);
	  for (int x = 0; x < X; ++x) {
		fprintf(fout, "%.2f,", map.Eout(x, y));
	  }
	  fprintf(fout, "\n");
	}
	*/


    // next section
    std::string oname3(argv[3]);
    oname3 += "stat_all_centres_" + timePeriodString + fileFormat;
    FILE *fout3 = fopen(oname3.c_str(), "w");

    //fprintf(fout, "\n\n--------------------- %s ---------------------\n", "Stat: all centres (x, y)");
    fprintf(fout3, "время, с||расстояние от центра, клеток\n");
    // { plot fix
    //fprintf(fout3, "0,%d,%d\n", 0, 0);
    //fprintf(fout3, "0,%d,%d\n", 0, Y);
    //fprintf(fout3, "0,%d,%d\n", X, 0);
    //fprintf(fout3, "0,%d,%d\n", X, Y);
    // } fix

    // (t, R_mean, R_sigma)
    std::list<std::tuple<float, float, float, float>> group_centres;
    std::list<float> current_Rs;

    float current_t = std::get<0>(all_centres_repeat.front());

    double xLC, yLC, distanceFromCenter;
    double xCenter = 100., yCenter = 100.;

    for (auto &&c : all_centres_repeat) {

        xLC = std::get<1>(c);
        yLC = std::get<2>(c);

        distanceFromCenter = sqrt((xLC - xCenter) * (xLC - xCenter) + (yLC - yCenter) * (yLC - yCenter));

        fprintf(fout3, "%.3f,%.3f\n", std::get<0>(c), distanceFromCenter);
        //fprintf(fout3, "%.3f, %.3f, %.3f\n", std::get<0>(c), xLC, yLC);

        
        if (std::get<0>(c) - current_t > 0.1f) { // calculate and drop
            auto M_D_SD = Get_M_D_SD(current_Rs);
            group_centres.emplace_back(current_t, std::get<0>(M_D_SD), std::get<1>(M_D_SD), std::get<2>(M_D_SD));

            current_t = std::get<0>(c);
            current_Rs.clear();
        }
        current_Rs.emplace_back(r_points(Point(100, 100), Point(std::get<1>(c), std::get<2>(c))));
    }
    //fclose(fout3);

    // next section
    std::string oname4(argv[3]);
    oname4 += "Rdispersion_" + timePeriodString + fileFormat;
    FILE *fout4 = fopen(oname4.c_str(), "w");
    fprintf(fout4, "время, с||среднее расстояние от центра, клеток||дисперсия, клеток||среднеквадратичное отклонение, клеток\n");
    for (auto &&c : group_centres) {
        fprintf(fout4, "%.3f,%.3f,%.3f,%.3f\n",
        std::get<0>(c), std::get<1>(c), std::get<2>(c), std::get<3>(c));
    }

    fclose(fout);
    fclose(fout2);
    fclose(fout3);
    fclose(fout4);

#endif
    std::cout << "Num cycles = " << numCycles << ", first ignition, sec = " <<  first_T << ", last ignition, sec = " << last_T << "\n";
    return 0;
}
