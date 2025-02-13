<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">
<html>
<head>

</head>

<!-- one2newt  THEN column/index.php -->

<body text="silver" LINK="silver" VLINK="silver" BACKGROUND="../home/Home_files/backgnd4.jpg" style="background-attachment: fixed">
<font size="3" face="arial" color="white">

<!-- http://localhost/~sjc/column/index.php  TO VIEW LOCALLY -->

<!-- do in php as javascript a pain -->
<!-- freqs from  http://spec.jpl.nasa.gov/ftp/pub/catalog/c028001.cat -->
<!-- As from http://cdsarc.u-strasbg.fr/viz-bin/Cat?J/A%2bAS/117/557  - cms96  SEARCH VIA ROTATION QUANTUM NUMBERS AND SETTING Vibrational quantum number TO 0 -->
<?
    if (isset($_POST['mol'])) { // test for existence of form parameters
    $mol = $_POST['mol'];

    $T1 = $_POST['T1'];
    $z = $_POST['z'];
    $tau  = $_POST['tau'];
    $v = $_POST['v'];
    
    $J = $_POST['J'];
    $red = $_POST['red'];
   // $frest= $_POST['frest'];
   // $freq= $_POST['freq'];

 
    $Z = $z + 1;
    $T = ($red == "C") ? $T1 : ($Z * 2.7)+($T1 - 3); //something fishy - option R always reset
 

    $kT = $T*1.381e-23;

   // printf(" %f %e <br>\n", $T, $kT);
    $h = 6.626e-34;

    if ($mol==co) { 
      $E1 = 3.8450*1.9878e-23;
      $E = (0.5*$J*($J+1))*$E1;

      if ($J==0) { 
        $frest = 115.271e9; 
        $A = 7.7e-8;       
      }                     
	
     if ($J==1) { 
        $frest = 230.538e9;
        $A = 7.4e-7;
}
      if ($J==2) { 
        $frest = 345.7959899e9;
        $A = 2.7e-06;
}

	if ($J==3) { 
        $frest = 461.0407682e9;
        $A = 6.5e-06;
}

	if ($J==4) { 
        $frest = 576267.9305e6;
        $A = 1.3e-05;
}

       if ($J==5) { 
        $frest = 691473.0763e6;
        $A = 2.3e-05;
}
 	if ($J==6) { 
        $frest = 806651.8060e6;
        $A = 3.7e-05;
}

	if ($J==7) { 
        $frest = 921799.7000e6;
        $A =  5.5e-05;
}	

if ($J==8) { 
        $frest = 1036912.3930e6;
        $A =  7.8e-05;
}

if ($J==9) { 
        $frest = 1151985.4520e6;
        $A = 1.1e-04;
}

}
    
if ($mol==hcn) { 
      

      $E1 = 2.9564*1.9878e-23;
      $E = (0.5*$J*($J+1))*$E1;
     // $mu = 2.98; //Debyes 
	$mu = 2.942; //http://spec.jpl.nasa.gov/ftp/pub/catalog/doc/d027003.cat

	if ($J==0) { 
        $frest = 88.6318470e9; 
       $A = (1.165e-11)*(pow(($frest*1e-9),3))*(pow($mu,2))*($J+1)/(2*$J+3);
	}
	if ($J==1){
     	$frest = 177.2612230e9;
	$A = (1.165e-11)*(pow(($frest*1e-9),3))*(pow($mu,2))*($J+1)/(2*$J+3);
	}
	 if ($J==2){
    	$frest = 265.8861800e9;
	$A = (1.165e-11)*(pow(($frest*1e-9),3))*(pow($mu,2))*($J+1)/(2*$J+3);
	}
 	if ($J==3){
   	  $frest = 354.5054759e9;
	$A = (1.165e-11)*(pow(($frest*1e-9),3))*(pow($mu,2))*($J+1)/(2*$J+3);
	}
 	if ($J==4){
    	 $frest = 443.1161554e9;
	$A = (1.165e-11)*(pow(($frest*1e-9),3))*(pow($mu,2))*($J+1)/(2*$J+3);
	}
	 if ($J==5){
   	  $frest = 531.7163875e9;
	$A = (1.165e-11)*(pow(($frest*1e-9),3))*(pow($mu,2))*($J+1)/(2*$J+3);
	}
	 if ($J==6){
     	$frest = 620.30409525e9;
	$A = (1.165e-11)*(pow(($frest*1e-9),3))*(pow($mu,2))*($J+1)/(2*$J+3);
	}
	 if ($J==7){
    	 $frest = 708.8772081e9;
	$A = (1.165e-11)*(pow(($frest*1e-9),3))*(pow($mu,2))*($J+1)/(2*$J+3);
	}
 

	}

if ($mol==hco) { 

   $E1 = 2.9750*1.9878e-23;
   $E = (0.5*$J*($J+1))*$E1;
  
  //$mu = 4.48; // Debyes - see MGB's sheet and Rohlfs p379   
    $mu = 3.888; //http://spec.jpl.nasa.gov/ftp/pub/catalog/doc/d029002.cat	

	 if ($J==0){
 
   	 $frest = 89.188518e9; 
  	  $A = (1.165e-11)*(pow(($frest*1e-9),3))*(pow($mu,2))*($J+1)/(2*$J+3); 
	}
 	if ($J==1){
 	   $frest = 178.3750650e9; 
 	 $A = (1.165e-11)*(pow(($frest*1e-9),3))*(pow($mu,2))*($J+1)/(2*$J+3);
	}
 	 if ($J==2){
   $frest = 267.5576190e9; 
 	 $A = (1.165e-11)*(pow(($frest*1e-9),3))*(pow($mu,2))*($J+1)/(2*$J+3);
	}
  	if ($J==3){
  
  	 $frest = 356.7342880e9; 
 	 $A = (1.165e-11)*(pow(($frest*1e-9),3))*(pow($mu,2))*($J+1)/(2*$J+3);
	}
  	if ($J==4){
  
  	 $frest = 445.9029960e9; 
 	 $A = (1.165e-11)*(pow(($frest*1e-9),3))*(pow($mu,2))*($J+1)/(2*$J+3);
	}
 	 if ($J==5){
  
  	 $frest = 535.0617755e9; 
 	 $A = (1.165e-11)*(pow(($frest*1e-9),3))*(pow($mu,2))*($J+1)/(2*$J+3);
	}
  	if ($J==6){
  
  	 $frest = 624.2086733e9; 
 	 $A = (1.165e-11)*(pow(($frest*1e-9),3))*(pow($mu,2))*($J+1)/(2*$J+3);
	}
  	if ($J==7){
  
  	 $frest = 713.3420900e9; 
 	 $A = (1.165e-11)*(pow(($frest*1e-9),3))*(pow($mu,2))*($J+1)/(2*$J+3); // spontanoeus EMISSION coefficient 
	}
 


}

if ($mol==cs) { 

  $E1 = 1.6342*1.9878e-23;
  $E = (0.5*$J*($J+1))*$E1;

 	
	if ($J==0){
 	$frest = 48.9909780e9;
  	$A = 1.74e-06;
 	}
  	
 	if ($J==1){
    	$frest = 97.9809500e9;
  	$A = 1.68e-05;
 	}

 	if ($J==2){
   	$frest = 146.9690330e9;
  	$A = 6.05e-05;
 	}

  	if ($J==3){
   	$frest = 195.9542260e9;
  	$A = 1.49e-04;
 	}

  	if ($J==4){
   	$frest = 244.93564350e9;
  	$A = 2.97e-04;
 	}

  	if ($J==5){
   	$frest = 293.9122440e9;
  	$A = 5.21e-04;
 	}

  	if ($J==6){
   	$frest = 342.8830000e9;
  	$A = 8.37e-04;
 	}

     	if ($J==7){
    	$frest = 391.8470300e9;
   	$A = 1.26e-03;          
 	}

}

// printf("J = %d   A= %e  freq = %e E1 = %e   kT = %e <br> Freq not being reset properly<br> J appears to be okay<br> retaining last molecule in program\n", $J,  $A,  $frest,  $E1, $kT);

    $blah = exp($E/$kT)/(1-(exp((-1*$h*$frest)/($kT))));
    $gu = 2*($J+1)+1;
    $g = (2*$J)+1;

    $x = -1*($E1/$kT);
    $Q = 1;

    $j = 1;

    while (((2*$j)+1)*exp(0.5*$j*($j+1)*$x) > 0.001){ 
       
     // printf("%f  %e %e  %e <br>\n", $j,  $x,  $E1, $kT);
      $Q = $Q + (((2*$j)+1)*exp(0.5*$j*($j+1)*$x));
      $j = $j+1;
      //printf("%f  \n", ((2*$j)+1)*exp(0.5*$j*($j+1)*$x));
    }


$f = $Q*$blah;
       $c = 3e5;

$N = (8*3.1416*(pow($frest,3))*$f*(1.06*$tau*$v))/((pow($c,3))*$A*$gu*1e10);

} else {       //must go at bottom or form will submit itself!
    $z = $T1 = $tau = $v = 0; 
    $red = "C";
     $mol = "co";
	$J = "0";
   
   // $frest = $freq = $N = $temp = 0;
	$N = $temp = 0;	
}
$mol_select = array("co", "hcn", "hco", "cs");
  $mol_select[$mol] = "selected";
$J_select = array("0", "1", "2", "3", "4", "5", "6", "7"); //so pull downs remember values
  $J_select[$J] = "selected";
?>

<p>
<center><br>
<font size=5 face="arial" >Column Density Calculator</font><br>
<p>
<form method="post" action="<? $_SERVER['PHP_SELF'] ?>#start">
</center>

<font size=2 face="arial">
<strong>  

This program will calculate the total column density from the
rotational transitions (up to J = 7-->8) for the most commonly
observed interstellar molecules (see Section 3.1 of <b> <a
href="../home/MD1376-final.pdf"><i>Curran et al., 2004, MNRAS 352,563 </a></i></b>). <p>

<center>


<!--http://oopweb.com/HTML/Documents/Writing%20HTML/VolumeFrames.html-->
<font size=3 face="arial">
<table border=0 cellpadding=12 cellspacing=6>
<tr>

<td valign=top align=left><b>
<strong>Molecule </b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<select name="mol">
    <option value="co"  <?= $mol_select["co"]  ?>>CO</option>
    <option value="hcn" <?= $mol_select["hcn"] ?>>HCN</option>
    <option value="hco" <?= $mol_select["hco"] ?>>HCO<sup>+</sup></option>
    <option value="cs"  <?= $mol_select["cs"]  ?>>CS</option>
  </select>

</td>
<a name = start>
<td valign=top nowrap>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<strong>J = <select name="J">

  // try looping this next bit to take CO up to 9 --> 10
 // for ($i=0; $i<=5; $i++){
 // for (int i=0; i<=5; i++){
  
 <!--  <option value= i <?= $J_select[i] ?>>i --> i+1</option> -->
 <!-- <option value= $i <?= $J_select[$i] ?>>$i --> $(i+1)</option> -->  // need this to lose - HTML issue?
//  }  // JUST NOT WORKING - STUCK WITH ADVERTISING CO AS 7 -> 8

	<option value= "0" <?= $J_select["0"] ?>>0 --> 1</option>
	<option value= "1" <?= $J_select["1"] ?>>1 --> 2 </option>
	<option value= "2" <?= $J_select["2"] ?>>2 --> 3 </option>
	<option value= "3" <?= $J_select["3"] ?>>3 --> 4 </option>
	<option value= "4" <?= $J_select["4"] ?>>4 --> 5 </option>
	<option value= "5" <?= $J_select["5"] ?>>5 --> 6 </option>
	<option value= "6" <?= $J_select["6"] ?>>6 --> 7 </option>
	<option value= "7" <?= $J_select["7"] ?>>7 --> 8 </option>
   
</select>
</td>
</tr>


<tr><td valign=top nowrap>
<b>Optical depth
</b>
<input type="text" name="tau"  value = "<?= $tau ?>" size="5"> <strong>
</td>
<td valign=top nowrap>
<b>Line width
</b>
<input type="text" name="v" value = "<?= $v ?>" size="5"><strong> km/s
</td>
</tr>

<td nowrap><b>Excitation <br>temperature</b>&nbsp;
<input type="text" name="T1"  value = "<?= $T1 ?>" size="5"> <strong>K
</td>


<td nowrap>
<b>Redshift dependent</b>
<? if ($red == "C") { ?>
<input type="radio" name="red" value="C" checked> <b>no</b>
<input type="radio" name="red" value="R"> <b>yes</b>
 <? } else { ?>
<input type="radio" name="red" value="C" > <b>no</b>
<input type="radio" name="red" value="R" checked> <b>yes</b>
 <br>
<strong>
<font size=3 face="arial"><br>
Rest freq of <span id="frest">  <? printf("%.3f", $frest/1e9); ?></span><b>GHz <br>
 redshifted to <span id="freq"> <? printf("%.3f", $frest/($Z*1e9)); ?></span><b>GHz.
<? } ?>

</td>
</tr>
<tr>
<td></td>
<td valign=top nowrap>
<b>Redshift, z = </b><input type="text" name="z" value = "<?= $z ?>" size="5"> <strong>
</td>
</tr>
</table>


<input type="submit" value="Calculate"><p>
<strong>
<font size=4>
Total column density is, N =  <span id="N"> <? printf("%.4e", $N); ?></span>    <b>cm<sup>-2</sup>  and T<sub>ex</sub> = <span id="T"><? printf("%.1f", $T); ?> K</span></font><b><br>
<font size=2>Used Einstein A-coefficient of A = <span id="N"> <? printf("%.4e", $A); ?></span>    <b>s<sup>-1</sup>
</form>
<font size=3>
<!-- <strong>
<center><a href="oh.php" target="pers"><font color="gold">OH transitions</font></a></strong></center>-->
<font size=2>

<p>
<hline>
<center><a href="http://spec.jpl.nasa.gov/ftp/pub/catalog/catdir.html" target="TARGET=fav_96706010">NASA JPL Catalogue</a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="http://www.physics.nist.gov/cgi-bin/micro/table5/start.pl" target="TARGET=fav_96706010">
NIST Recommended Rest Frequencies for Observed Interstellar Molecular Microwave Transitions (LOVAS)</a><br>
<font size=2>

<p>
<font size=2>
<strong>  
If there are other molecules/transitions you'd like
included, just <a href="mailto:sjc@phys.unsw.edu.au?subject=Column
density calculator">send me </a>the Einstein A-coefficients/dipole moments.<b><br>


<strong>
<center><a href="oh.php" target="pers">OH transitions</font></a>  &nbsp;&nbsp;&nbsp;&nbsp <a href="../home/programs.html" target="pers">Back</a></strong>
<p></center>

</body>
</html>
