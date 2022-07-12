     #include <plotter.h>
     #include <fstream>
     const int maxorder = 12;
     
     void draw_c_curve (Plotter& plotter, double dx, double dy, int order)
     {
       if (order >= maxorder)
         plotter.fcontrel (dx, dy);	// continue path along (dx, dy)
       else
         {
           draw_c_curve (plotter,
                         0.5 * (dx - dy), 0.5 * (dx + dy), order + 1);
           draw_c_curve (plotter,
                         0.5 * (dx + dy), 0.5 * (dy - dx), order + 1);
         }
     }
     
     int main ()
     {
       // set a Plotter parameter
       PlotterParams params;
       params.setplparam ("PAGESIZE", (char *)"letter");
     
       std::ofstream outfile;
       outfile.open("plotout_c.ps",std::ios::out);
       PSPlotter plotter(cin, outfile, cerr, params); // declare Plotter
       if (plotter.openpl () < 0)                  // open Plotter
         {
           cerr << "Couldn't open Plotter\n";
           return 1;
         }
     
       plotter.fspace (0.0, 0.0, 1000.0, 1000.0); // specify user coor system
       plotter.flinewidth (0.25);       // line thickness in user coordinates
       plotter.pencolorname ("red");    // path will be drawn in red
       plotter.erase ();                // erase Plotter's graphics display
       plotter.fmove (600.0, 300.0);    // position the graphics cursor
       draw_c_curve (plotter, 0.0, 400.0, 0);
       if (plotter.closepl () < 0)      // close Plotter
         {
           cerr << "Couldn't close Plotter\n";
           return 1;
         }
       return 0;
     }
