#include <fstream>
#include <iostream>
using namespace std;

void tower_pos_gen()
{

  // The Default unit if the fun4all is cm
  // Here are all the parameters you would like to change
  // One thing need to be noticed here, when you want to change the tower width,
  // no matter the inner or outer tower, please keep their width ratio always at 2 (outer/inner)
  // Or there will be an overlap problem.
  // material para: 0.5[G4_PbWO4],
  //                1.5[G4_GLASS_LEAD],
  //                2.5[G4_BARIUM_SULFATE],
  //                3.5[G4_CESIUM_IODIDE],
  //                4.5[Sci-glass-selfdefined]
  //                5.5[heavier Sci-glass-selfdefined]
  //                6.5[Sci-glass(lead)-selfdefined]
  //                7.5[Sci-glass Nathaly]
  //                8.5[Sci-glass g4e]
  
  double PWO_carbon_Gap = 0.1,
         PWO_air_Gap = 0.,
         PWO_Gap = PWO_carbon_Gap + PWO_air_Gap,
         PWO_tower_width = 2.0,
         PWO_Width = (PWO_tower_width - 2. * PWO_Gap),
         PWO_Thickness = 20,
         PWO_OuterR = 53.,
         PWO_InnerR = 10.,
         Inner_mat = 0.5; 
  
  double Glass_carbon_Gap = 0.1,
         Glass_air_Gap = 0.,
         Glass_Gap = Glass_carbon_Gap + Glass_air_Gap,
         Glass_tower_width = PWO_tower_width * 2.,
         Glass_Width = (Glass_tower_width - 2. * Glass_Gap),
         Glass_Thickness = 40.,
         Glass_InnerR = PWO_OuterR,
         Glass_OuterR = 85.,
         Outer_mat = 7.5; 
  
  double inner_glo_zpos = -200.,
         outer_glo_zpos = inner_glo_zpos - (Glass_Thickness - PWO_Thickness) / 2.;

  char Y_N;


  /*
  cout << endl << " The following are default parameters for the EEMC configuartion......." << endl << endl;
  cout << "=================================================================" << endl << endl;
  cout << " Inner configuration:" << endl;
  cout << " inner R[cm]: " << PWO_InnerR << "  ||  outer R[cm]: " << PWO_OuterR << endl;
  cout << " crystal width[cm]: " << PWO_Width;
  cout << "  ||  shell thickness[cm]: " << PWO_carbon_Gap;
  cout << "  ||  air gap[cm]: " << PWO_air_Gap << endl;
  cout << " crystal thickness[cm]: " << PWO_Thickness << endl;
  cout << " Global position of inner: [0, 0, " << inner_glo_zpos << "cm]" << endl << endl;
  cout << " Outer configuration:" << endl;
  cout << " inner R[cm]: " << Glass_InnerR << "  ||  outer R[cm]: " << Glass_OuterR << endl;
  cout << " sci-glass width[cm]: " << Glass_Width;
  cout << "  ||  shell thickness[cm]: " << Glass_carbon_Gap;
  cout << "  ||  air gap[cm]: " << Glass_air_Gap << endl;
  cout << " sci-glass thickness[cm]: " << Glass_Thickness << endl;
  cout << " Global position of outer: [0, 0, " << outer_glo_zpos << "cm]" << endl << endl;
  cout << "=================================================================" << endl << endl << endl;

  cout << " If you want to personalize the ones for yourself, please press [Y], or press [N] for default value......" << endl << endl;
  cin >> Y_N;
  cout << endl << endl;
 
  if( (Y_N == 'Y') || (Y_N == 'y') )
    {
      cout << " Please key in the desired parameters for the inner configuration......" << endl;
      cout << " inner R[cm]: ";
      cin >> PWO_InnerR;
      cout << endl << " outer R[cm]: ";
      cin >> PWO_OuterR;
      cout << endl << " crystal width[cm]: ";
      cin >> PWO_Width;
      cout << endl << " shell thickness[cm]: ";
      cin >> PWO_carbon_Gap;
      cout << endl << " air gap[cm]: ";
      cin >> PWO_air_Gap;
      cout << "Tower width for the inner is: " << PWO_Width + 2. * (PWO_carbon_Gap + PWO_air_Gap) << endl;
      cout << "Global z position of inner, should be < 0: ";
      cin >> inner_glo_zpos;
      cout << endl << endl;
      
      cout << " Please key in the desired parameters for the outer configuration......" << endl;
      Glass_InnerR = PWO_OuterR ;
      cout << endl << " outer R[cm]: ";
      cin >> Glass_OuterR;
      cout << endl << " shell thickness[cm]: ";
      cin >> Glass_carbon_Gap;
      cout << endl << " air gap[cm]: ";
      cin >> Glass_air_Gap;
      Glass_Width = (PWO_Width + 2. * (PWO_carbon_Gap + PWO_air_Gap)) * 2. - 2. * (Glass_carbon_Gap + Glass_air_Gap);
      cout << "Tower width for the inner is: " << Glass_Width + 2. * (Glass_carbon_Gap + Glass_air_Gap) << endl;
      cout << "Global z position of inner, should be < 0: ";
      cin >> outer_glo_zpos;
    }
  else if( (Y_N == 'N') || (Y_N == 'n') )
    cout << "The default parameters would be applied...." << endl;
  else
    cout << "Meaning less choice, the default value would be apply." << endl << endl;

  //         outer_glo_zpos = inner_glo_zpos + (PWO_Thickness / 2. - Glass_Thickness / 2.);
  */


  

  
  // CRYSTAL
  double diameter = 2 * PWO_OuterR;
  double leftTowerPos, topTowerPos;
  // How many towers do we have per row/column?
  // Add a gap + diameter as if we have N towers, we have N-1 gaps;
  //  int towersInRow = std::ceil((diameter + 2. * PWO_Gap) /  (PWO_Width + 2. * PWO_Gap));
  //  int towersInRow = std::floor((diameter + 2. * PWO_Gap) /  (PWO_Width + 2. * PWO_Gap));
  int towersInRow = std::floor((diameter + 2. * Glass_Gap) /  (Glass_Width + 2. * Glass_Gap));
  //  cout << "towersInRow: " << towersInRow << endl;
  
  // Is it odd or even number of towersInRow

  if(towersInRow%2)
    {
      //             |
      //      [ ][ ][ ][ ][ ]
      //       ^     |
      int towersInHalfRow = std::ceil(towersInRow/2.0);
      //      leftTowerPos = topTowerPos = -towersInHalfRow * (PWO_Width + 2. * PWO_Gap);
      leftTowerPos = topTowerPos = -towersInHalfRow * (Glass_Width + 2. * Glass_Gap);
    }
  else
    {
      //               |
      //      [ ][ ][ ][ ][ ][ ]
      //       ^      |
      int towersInHalfRow = towersInRow/2;
      //      leftTowerPos = topTowerPos = -(towersInHalfRow - 0.5) * (PWO_Width + 2. * PWO_Gap);
      leftTowerPos = topTowerPos = -(towersInHalfRow - 0.5) * (Glass_Width + 2. * Glass_Gap);
    }

  std::printf("\nCE EMCAL PWO SQUARE START\n");
  std::printf("PWO_Thickness  = %.1f cm;\n", PWO_Thickness );
  std::printf("PWO_Width      = %.1f cm;\n", PWO_Width );
  std::printf("PWO_total_Gap  = %.1f cm;\n", PWO_Gap );
  std::printf("PWO_carbon_Gap = %.1f cm;\n", PWO_carbon_Gap );
  std::printf("PWO_air_Gap    = %.1f cm;\n", PWO_air_Gap );
  std::printf("PWO_InnerR     = %.1f cm;\n", PWO_InnerR );
  std::printf("PWO_OuterR     = %.1f cm;\n", PWO_OuterR );
  //  std::printf("PWO_PosZ      = %.1f cm;\n", PWO_PosZ );
  std::printf("Towers in Row/Col   = %d;\n", towersInRow);
  std::printf("Top left tower pos  = %.1f, %.1f cm;\n", leftTowerPos , topTowerPos );


  int towerIndex = 0, count = 0;
  ofstream ofs;

  
  // =========================================
  // Generating the crystal position mapping
  // =========================================
  
  ofs.open("tower_map_crystal.txt");
  if (!ofs.is_open()) 
    cout << "Failed to open file." << endl;

  //The prefix parameter of the mapping
  ofs << "Gtype 1" << endl;
  ofs << "Gr1_inner " << PWO_InnerR << endl;
  ofs << "Gr1_outer " << PWO_OuterR << endl;
  ofs << "Gr2_inner " << PWO_InnerR << endl;
  ofs << "Gr2_outer " << PWO_OuterR << endl;
  ofs << "Gdz " << PWO_Thickness << endl;
  ofs << "Gx0 " << 0 << endl;
  ofs << "Gy0 " << 0 << endl;
  ofs << "Gz0 " << 0 << endl;
  ofs << "Grot_x " << 0 << endl;
  ofs << "Grot_y " << 0 << endl;
  ofs << "Grot_z " << 0 << endl;
  ofs << "Gcrystal_dx " << PWO_Width << endl;
  ofs << "Gcrystal_dy " << PWO_Width << endl;
  ofs << "Gcrystal_dz " << PWO_Thickness << endl;
  ofs << "Gcarbon_gap " << PWO_carbon_Gap << endl;
  ofs << "Gair_gap " << PWO_air_Gap << endl;
  ofs << "Gmaterial " << Inner_mat << endl;
  ofs << "Gcolor_R " << 0.0 << endl;
  ofs << "Gcolor_G " << 0.0 << endl;
  ofs << "Gcolor_B " << 1.0 << endl;
  ofs << "#Tower: type, idx_j, idx_k, idx_l, x[cm], y[cm], z[cm], dx[cm], dy[cm], dz[cm], rot_x, rot_y, rot_z" << endl;
  
  for(int colIndex = 0 ; colIndex < towersInRow ; colIndex++)
    {
      for(int rowIndex = 0 ; rowIndex < towersInRow ; rowIndex++)
	{
	  // double x = leftTowerPos + colIndex * (PWO_Width + 2. * PWO_Gap);
	  // double y = topTowerPos + rowIndex * (PWO_Width + 2. * PWO_Gap);
	  double x = leftTowerPos + colIndex * (Glass_Width + 2. * Glass_Gap);
	  double y = topTowerPos + rowIndex * (Glass_Width + 2. * Glass_Gap);
	  double z = inner_glo_zpos;
	  //	  double z = -250.;
	  double Cry_R = std::sqrt(x * x + y * y);
	  double inner_r = 0.;
	  
	  //	  if ( (std::abs(y) < PWO_OuterR && std::abs(x) < PWO_OuterR) && (std::abs(y) > PWO_InnerR || std::abs(x) > PWO_InnerR) )
	  //	  if ( (Cry_R < PWO_OuterR) && (std::abs(y) > (PWO_InnerR + 2.) || std::abs(x) > (PWO_InnerR + 2.) ) )
	  if ( (Cry_R < PWO_OuterR) && (Cry_R > PWO_InnerR) )
	    {
	      int index_dis = 0;
	      
	      for(int i = 0 ; i < 2 ; i++)
		{
		  for(int j = 0 ; j < 2 ; j++)
		    {
		      double x_cry, y_cry;
		      int cry_col_index, cry_row_index;

		      if( (i==0) && (j==0) )
			{
			  x_cry = x - (PWO_Width + 2. * PWO_Gap) / 2.;
			  y_cry = y - (PWO_Width + 2. * PWO_Gap) / 2.;
			  inner_r = std::sqrt(x_cry * x_cry + y_cry * y_cry);
			  if(inner_r < PWO_InnerR)
			    continue;
			  cry_col_index = colIndex * 2 + i;
			  cry_row_index = rowIndex * 2 + j;
			}
		      else if( (i==0) && (j==1) )
			{
			  x_cry = x - (PWO_Width + 2. * PWO_Gap) / 2.;
			  y_cry = y + (PWO_Width + 2. * PWO_Gap) / 2.;
			  inner_r = std::sqrt(x_cry * x_cry + y_cry * y_cry);
			  if(inner_r < PWO_InnerR)
			    continue;
			  cry_col_index = colIndex * 2 + i;
			  cry_row_index = rowIndex * 2 + j;
			}
		      else if( (i==1) && (j==0) )
			{
			  x_cry = x + (PWO_Width + 2. * PWO_Gap) / 2.;
			  y_cry = y - (PWO_Width + 2. * PWO_Gap) / 2.;
			  inner_r = std::sqrt(x_cry * x_cry + y_cry * y_cry);
			  if(inner_r < PWO_InnerR)
			    continue;
			  cry_col_index = colIndex * 2 + i;
			  cry_row_index = rowIndex * 2 + j;
			}
		      else if( (i==1) && (j==1) )
			{
			  x_cry = x + (PWO_Width + 2. * PWO_Gap) / 2.;
			  y_cry = y + (PWO_Width + 2. * PWO_Gap) / 2.;
			  inner_r = std::sqrt(x_cry * x_cry + y_cry * y_cry);
			  if(inner_r < PWO_InnerR)
			    continue;
			  cry_col_index = colIndex * 2 + i;
			  cry_row_index = rowIndex * 2 + j;
			}

		      int code = 1000 * cry_col_index + cry_row_index;

		      ofs << "Tower " << 0 << " ";
		      ofs << cry_col_index << " " << cry_row_index << " " << 0 << " ";
		      ofs << x_cry << " " << y_cry << " " << z << " ";
		      ofs << PWO_Width << " " << PWO_Width << " " << PWO_Thickness << " ";
		      ofs << 0 << " " << 0 << " " << 0 << endl;
	      
		      towerIndex++;
		      count++;
		      
		    }
		} // replace the glass pos to the crystal ones
	    } // inner range 
	  
	}
    }

  ofs.close();
  towerIndex = 0;
  cout << "Inner use : " << count << " towers to construct!!!" << endl << endl;
  count = 0;








 
  // =========================================
  // Generating the Glass position mapping
  // =========================================

  // Glass
  diameter = 2 * Glass_OuterR;
  // How many towers do we have per row/column?
  // Add a gap + diameter as if we have N towers, we have N-1 gaps;
  //  towersInRow = std::ceil((diameter + 2. * Glass_Gap) /  (Glass_Width + 2. * Glass_Gap));
  towersInRow = std::floor((diameter + 2. * Glass_Gap) /  (Glass_Width + 2. * Glass_Gap));

  // Is it odd or even number of towersInRow

  if(towersInRow%2)
    {
      //             |
      //      [ ][ ][ ][ ][ ]
      //       ^     |
      int towersInHalfRow = std::ceil(towersInRow/2.0);
      leftTowerPos = topTowerPos = -towersInHalfRow * (Glass_Width + 2. * Glass_Gap);
    }
  else
    {
      //               |
      //      [ ][ ][ ][ ][ ][ ]
      //       ^      |
      int towersInHalfRow = towersInRow/2;
      leftTowerPos = topTowerPos = -(towersInHalfRow - 0.5) * (Glass_Width + 2. * Glass_Gap);
    }

  std::printf("\nCE EMCAL Glass SQUARE START\n");
  std::printf("Glass_Thickness   = %.1f cm;\n", Glass_Thickness );
  std::printf("Glass_Width       = %.1f cm;\n", Glass_Width );
  std::printf("Glass_Gap         = %.1f cm;\n", Glass_Gap );
  std::printf("Glass_carbon_Gap  = %.1f cm;\n", Glass_carbon_Gap );
  std::printf("Glass_air_Gap     = %.1f cm;\n", Glass_air_Gap );
  std::printf("Glass_InnerR      = %.1f cm;\n", Glass_InnerR );
  std::printf("Glass_OuterR      = %.1f cm;\n", Glass_OuterR );
  //  std::printf("Glass_PosZ      = %.1f cm;\n", Glass_PosZ );
  std::printf("Towers in Row/Col   = %d;\n", towersInRow);

  std::printf("Top left tower pos  = %.1f, %.1f cm;\n", leftTowerPos , topTowerPos );

  
  ofs.open("tower_map_glass.txt");
  if (!ofs.is_open()) 
    cout << "Failed to open file." << endl;
  
  ofs << "Gtype 1" << endl;
  ofs << "Gr1_inner " << Glass_InnerR << endl;
  ofs << "Gr1_outer " << Glass_OuterR << endl;
  ofs << "Gr2_inner " << Glass_InnerR << endl;
  ofs << "Gr2_outer " << Glass_OuterR << endl;
  ofs << "Gdz " << Glass_Thickness << endl;
  ofs << "Gx0 " << 0 << endl;
  ofs << "Gy0 " << 0 << endl;
  ofs << "Gz0 " << 0 << endl;
  ofs << "Grot_x " << 0 << endl;
  ofs << "Grot_y " << 0 << endl;
  ofs << "Grot_z " << 0 << endl;
  ofs << "Gcrystal_dx " << Glass_Width << endl;
  ofs << "Gcrystal_dy " << Glass_Width << endl;
  ofs << "Gcrystal_dz " << Glass_Thickness << endl;
  ofs << "Gcarbon_gap " << Glass_carbon_Gap << endl;
  ofs << "Gair_gap " << Glass_air_Gap << endl;
  ofs << "Gmaterial " << Outer_mat << endl;
  ofs << "Gcolor_R " << 1.0 << endl;
  ofs << "Gcolor_G " << 0.0 << endl;
  ofs << "Gcolor_B " << 0.0 << endl;
  ofs << "#Tower: type, idx_j, idx_k, idx_l, x[cm], y[cm], z[cm], dx[cm], dy[cm], dz[cm], rot_x, rot_y, rot_z" << endl;
  
  for(int colIndex = 0 ; colIndex < towersInRow ; colIndex++)
    {
      for(int rowIndex = 0 ; rowIndex < towersInRow ; rowIndex++)
	{
	  double x = leftTowerPos + colIndex * (Glass_Width + 2. * Glass_Gap);
	  double y = topTowerPos + rowIndex * (Glass_Width + 2. * Glass_Gap);
	  double z = outer_glo_zpos;
	  //	  double z = -260.;
	  double Cry_R = std::sqrt(x * x + y * y);

	  //	  if ( (Cry_R < Glass_OuterR) && (std::abs(y) > Glass_InnerR || std::abs(x) > Glass_InnerR) )
	  if ( (Cry_R < Glass_OuterR) && (Cry_R > Glass_InnerR) ) // To avoid overlapped with its mother volume
	    {
	      int code = 1000000 + 1000 * rowIndex + colIndex;

	      ofs << "Tower " << 0 << " ";
	      ofs << colIndex << " " << rowIndex << " " << 0 << " ";
	      ofs << x << " " << y << " " << z << " ";
	      ofs << Glass_Width << " " << Glass_Width << " " << Glass_Thickness << " ";
	      ofs << 0 << " " << 0 << " " << 0 << endl;
	      
	      towerIndex++;
	      count++;
	    }
	}
    }

  cout << "Outer use : " << count << " towers to construct!!!" << endl << endl;
  ofs.close();

  
    
  //  cout << "Crystal construction result: ideal constructed: " << (towersInRow * towersInRow) << " || real constructed: " << towerIndex << endl;
      
  
}
