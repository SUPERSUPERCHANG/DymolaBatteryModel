within ;
package ETBATTERY
  model ParameterLookup
    Modelica.Blocks.Tables.CombiTable2Ds R1C(
      tableOnFile=true,
      table=[0.0,-20,0.0; 0.1345,0.11433987088379959,0.0; 0.0,0.0,0.0],
      tableName="R12_SOC_c",
      fileName="D:/Battery/ET battery/matlab.mat",
      smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
      extrapolation=Modelica.Blocks.Types.Extrapolation.NoExtrapolation)
      annotation (Placement(transformation(extent={{-16,50},{4,70}})));
    Modelica.Blocks.Tables.CombiTable2Ds R1D(
      tableOnFile=true,
      table=[0.0,0.0; 0.0,0.0],
      tableName="R12_SOC_d",
      fileName="D:/Battery/ET battery/matlab.mat",
      smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
      extrapolation=Modelica.Blocks.Types.Extrapolation.NoExtrapolation)
      annotation (Placement(transformation(extent={{-16,16},{4,36}})));
    Modelica.Blocks.Interfaces.RealInput SOC annotation (Placement(
          transformation(extent={{-114,26},{-86,54}}), iconTransformation(
            extent={{-114,26},{-86,54}})));
    Modelica.Blocks.Interfaces.RealInput Ts annotation (Placement(
          transformation(extent={{-114,-14},{-86,14}}), iconTransformation(
            extent={{-114,-14},{-86,14}})));
    Modelica.Blocks.Interfaces.RealInput I annotation (Placement(transformation(
            extent={{-114,-54},{-86,-26}}), iconTransformation(extent={{-114,
              -54},{-86,-26}})));
    Modelica.Blocks.Logical.Switch switch1
      annotation (Placement(transformation(extent={{42,6},{62,26}})));
    Modelica.Blocks.Logical.GreaterEqualThreshold greaterEqualThreshold(
        threshold=0)
      annotation (Placement(transformation(extent={{-74,-26},{-62,-14}})));
    Modelica.Blocks.Interfaces.RealOutput R1 annotation (Placement(
          transformation(extent={{88,16},{112,40}}), iconTransformation(extent=
              {{88,16},{112,40}})));
    Modelica.Blocks.Tables.CombiTable2Ds R2C(
      tableOnFile=true,
      table=[0.0,0.0,0.0],
      tableName="R22_SOC_c",
      fileName="D:/Battery/ET battery/matlab.mat")
      annotation (Placement(transformation(extent={{-16,-40},{4,-20}})));
    Modelica.Blocks.Tables.CombiTable2Ds R2D(
      tableOnFile=true,
      table=[0.0,0.0; 0.0,0.0],
      tableName="R22_SOC_d",
      fileName="D:/Battery/ET battery/matlab.mat")
      annotation (Placement(transformation(extent={{-16,-80},{4,-60}})));
    Modelica.Blocks.Logical.Switch switch2
      annotation (Placement(transformation(extent={{42,-54},{62,-34}})));
    Modelica.Blocks.Interfaces.RealOutput R2 annotation (Placement(
          transformation(extent={{88,-52},{112,-28}}), iconTransformation(
            extent={{88,-52},{112,-28}})));
  equation
    connect(SOC, R1C.u1) annotation (Line(points={{-100,40},{-28,40},{-28,66},{
            -18,66}}, color={0,0,127}));
    connect(Ts, R1C.u2) annotation (Line(points={{-100,0},{-26,0},{-26,54},{-18,
            54}}, color={0,0,127}));
    connect(SOC, R1D.u1) annotation (Line(points={{-100,40},{-46,40},{-46,32},{
            -18,32}}, color={0,0,127}));
    connect(Ts, R1D.u2) annotation (Line(points={{-100,1.77636e-15},{-74,
            1.77636e-15},{-74,20},{-18,20}}, color={0,0,127}));
    connect(R1C.y, switch1.u1) annotation (Line(points={{5,60},{34,60},{34,24},
            {40,24}}, color={0,0,127}));
    connect(R1D.y, switch1.u3) annotation (Line(points={{5,26},{12,26},{12,0},{
            40,0},{40,8}}, color={0,0,127}));
    connect(I, greaterEqualThreshold.u) annotation (Line(points={{-100,-40},{
            -80,-40},{-80,-20},{-75.2,-20}}, color={0,0,127}));
    connect(greaterEqualThreshold.y, switch1.u2) annotation (Line(points={{
            -61.4,-20},{-56,-20},{-56,6},{32,6},{32,16},{40,16}}, color={255,0,
            255}));
    connect(switch1.y, R1) annotation (Line(points={{63,16},{78,16},{78,28},{
            100,28}}, color={0,0,127}));
    connect(R1, R1)
      annotation (Line(points={{100,28},{100,28}}, color={0,0,127}));
    connect(greaterEqualThreshold.y, switch2.u2) annotation (Line(points={{
            -61.4,-20},{-40,-20},{-40,-44},{40,-44}}, color={255,0,255}));
    connect(SOC, R2C.u1) annotation (Line(points={{-100,40},{-46,40},{-46,32},{
            -24,32},{-24,-24},{-18,-24}}, color={0,0,127}));
    connect(SOC, R2D.u1) annotation (Line(points={{-100,40},{-46,40},{-46,32},{
            -24,32},{-24,-64},{-18,-64}}, color={0,0,127}));
    connect(Ts, R2C.u2) annotation (Line(points={{-100,0},{-26,0},{-26,-36},{
            -18,-36}}, color={0,0,127}));
    connect(Ts, R2D.u2) annotation (Line(points={{-100,0},{-26,0},{-26,-76},{
            -18,-76}}, color={0,0,127}));
    connect(switch2.y, R2) annotation (Line(points={{63,-44},{82,-44},{82,-40},
            {100,-40}}, color={0,0,127}));
    connect(R2C.y, switch2.u1)
      annotation (Line(points={{5,-30},{40,-30},{40,-36}}, color={0,0,127}));
    connect(R2D.y, switch2.u3) annotation (Line(points={{5,-70},{34,-70},{34,
            -52},{40,-52}}, color={0,0,127}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(extent={{-54,66},{58,-68}}, lineColor={28,108,200}),
          Line(points={{-100,44}}, color={28,108,200}),
          Text(
            extent={{-52,40},{-34,40}},
            textColor={28,108,200},
            fontSize=72,
            textString="SOC"),
          Text(
            extent={{-48,50},{-16,22}},
            textColor={28,108,200},
            textString="SOC"),
          Text(
            extent={{-50,10},{-20,-10}},
            textColor={28,108,200},
            textString="Ts"),
          Text(
            extent={{-44,-34},{-26,-50}},
            textColor={28,108,200},
            textString="I"),
          Text(
            extent={{26,38},{56,18}},
            textColor={28,108,200},
            textString="R1"),
          Line(points={{-86,40},{-54,40}}, color={28,108,200}),
          Line(points={{-86,0},{-54,0}}, color={28,108,200}),
          Line(points={{-84,-40}}, color={28,108,200}),
          Line(points={{-86,-40},{-54,-40}}, color={28,108,200}),
          Line(points={{88,28},{58,28}}, color={28,108,200}),
          Text(
            extent={{26,-30},{52,-48}},
            textColor={28,108,200},
            textString="R2"),
          Line(points={{58,-40},{88,-40}}, color={28,108,200})}),  Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end ParameterLookup;

  model CCCV_Charging
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{-10,-98},{10,-78}})));
    Modelica.Electrical.Analog.Sensors.PowerSensor powerSensor
      annotation (Placement(transformation(extent={{-10,46},{10,66}})));
    Modelica.Blocks.Continuous.Integrator energy(u(unit="W"), y(unit="J"))
      annotation (Placement(transformation(extent={{38,8},{58,28}})));
    Modelica.Electrical.Batteries.Utilities.CCCVcharger
                          cccvCharger(I=50, Vend=42) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-66,6})));
    replaceable CylindricBat cylindricBat annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={18,-8})));
  equation
    connect(powerSensor.nv,ground. p)
      annotation (Line(points={{0,46},{0,-78}},     color={0,0,255}));
    connect(powerSensor.pc,powerSensor. pv) annotation (Line(points={{-10,56},{
            -16,56},{-16,70},{0,70},{0,66}},
                                color={0,0,255}));
    connect(powerSensor.power,energy. u)
      annotation (Line(points={{-10,45},{-10,38},{30,38},{30,18},{36,18}},
                                                           color={0,0,127}));
    connect(powerSensor.pc,cccvCharger. p) annotation (Line(points={{-10,56},{
            -66,56},{-66,16}},  color={0,0,255}));
    connect(ground.p,cccvCharger. n) annotation (Line(points={{0,-78},{0,-10},{
            -66,-10},{-66,-4}},
                             color={0,0,255}));
    connect(powerSensor.nc, cylindricBat.p)
      annotation (Line(points={{10,56},{18,56},{18,2}}, color={0,0,255}));
    connect(cylindricBat.n, ground.p) annotation (Line(points={{18,-18},{18,-72},
            {0,-72},{0,-78}}, color={0,0,255}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end CCCV_Charging;

  model CellRCStack
    "Battery with open-circuit voltage dependent on state of charge, self-discharge, inner resistance and a series of RC-elements"
    extends Modelica.Electrical.Batteries.BaseClasses.BaseCellStack(r0(final R=
            Ns*ParameterLookup.R12_SOC_C/Np),
      redeclare Modelica.Electrical.Batteries.ParameterRecords.TransientData.CellData cellData);
    extends Modelica.Electrical.Batteries.Icons.TransientModel;
    Modelica.Electrical.Analog.Basic.Resistor resistor[cellData.nRC](
      final R=Ns*cellData.rcData.R/Np,
      final T_ref=cellData.rcData.T_ref,
      final alpha=cellData.rcData.alpha,
      each final useHeatPort=true)
      annotation (Placement(transformation(extent={{30,-30},{50,-10}})));
    Modelica.Electrical.Analog.Basic.Capacitor capacitor[cellData.nRC](each v(
          fixed=true, each start=0), final C=Np*cellData.rcData.C/Ns)
      annotation (Placement(transformation(extent={{30,30},{50,10}})));
    ParameterLookup parameterLookup
      annotation (Placement(transformation(extent={{-90,-72},{-70,-52}})));
    ParameterLookup parameterLookup1 annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-114,60})));
    Modelica.Blocks.Interaction.Show.RealValue realValue annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-118,90})));
    Modelica.Blocks.Interaction.Show.RealValue realValue1 annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-110,90})));
  equation
    assert(cellData.R0 > 0, "Ri has to be greater than sum(rcParameters.R)");
    assert(cellData.rcData[1].R > 0, "Parameters of RC-elements undefined!");
    //connect the RC-elements
    connect(resistor[1].p, r0.n)
      annotation (Line(points={{30,-20},{30,0},{10,0}},color={0,0,255}));
    for k in 1:cellData.nRC loop
      connect(capacitor[k].p, resistor[k].p)
        annotation (Line(points={{30,20},{30,-20}}, color={0,0,255}));
      connect(capacitor[k].n, resistor[k].n)
        annotation (Line(points={{50,20},{50,-20}},          color={0,0,255}));
      connect(internalHeatPort, resistor[k].heatPort)
        annotation (Line(points={{0,-80},{0,-40},{40,-40},{40,-30}}, color={191,0,0}));
      if k < cellData.nRC then
        connect(resistor[k].n, resistor[k + 1].p);
      end if;
    end for;
    connect(resistor[cellData.nRC].n, n)
      annotation (Line(points={{50,-20},{50,0},{100,0}},color={0,0,255}));
    connect(limIntegrator.u, parameterLookup1.I) annotation (Line(points={{-80,18},
            {-80,14},{-110,14},{-110,50}},     color={0,0,127}));
    connect(ocv_soc.u, parameterLookup1.SOC) annotation (Line(points={{-72,50},
            {-80,50},{-80,46},{-82,46},{-82,44},{-104,44},{-104,42},{-108,42},{
            -108,40},{-120,40},{-120,50},{-118,50}},           color={0,0,127}));
    connect(parameterLookup1.R1, realValue.numberPort) annotation (Line(points={{-116.8,
            70},{-118,70},{-118,78.5}},                      color={0,0,127}));
    connect(parameterLookup1.R2, realValue1.numberPort)
      annotation (Line(points={{-110,70},{-110,78.5}}, color={0,0,127}));
    annotation (
      Documentation(info="<html>
<p>
Extends the model <a href=\"modelica://Modelica.Electrical.Batteries.BatteryStacks.CellStack\">CellStack</a> by a series of RC-elements, describing the transient behaviour of the battery.
</p>
<p>
This model can be used for a single cell <code>Ns = Np = 1</code> as well as a stack built from identical cells.
</p>
<p>
For details, see <a href=\"modelica://Modelica.Electrical.Batteries.UsersGuide.Concept\">concept</a> and <a href=\"modelica://Modelica.Electrical.Batteries.UsersGuide.Parameterization\">parameterization</a>.
</p>
<h4>Note</h4>
<p>
Parameter record array <a href=\"modelica://Modelica.Electrical.Batteries.ParameterRecords.TransientData.RCData\">rcData</a> contained in
parameter record <a href=\"modelica://Modelica.Electrical.Batteries.ParameterRecords.TransientData.CellData\">cellData</a> has to be specified.
</p>
<p>
The total inner resistance is the sum of the resistance of resistor <code>r0</code> and the sum of the resistances of the resistors of the RC-elements.
</p>
</html>"));
  end CellRCStack;

  model CylindricBat "Cell model template"

    extends Battery.Cells.Partials.PartialCylindricCell;
    extends Battery.Common.Icons.Boxed_BottomName;

    replaceable Battery.Cells.Thermal.Variants.CylindricConstTemp thermalModel constrainedby
      Battery.Cells.Thermal.Partials.PartialCylindric(
      final N_surface=N_surface, final N_verticalElements=N_verticalElements, final T_init=T_init) "Thermal model in cell"
      annotation (Placement(transformation(extent={{-20,40},{20,80}})), choicesAllMatching=true,
        Dialog(group="Replaceable Models"));

    replaceable Battery.Cells.Electric.Partials.PartialElectric electricModel constrainedby
      Battery.Cells.Electric.Partials.PartialElectric(final SOC_init=SOC_init, final T_init=T_init) "Electric model in cell"
      annotation (
      Dialog(group="Replaceable Models"),
      choicesAllMatching=true,
      Placement(transformation(extent={{-20,-20},{20,20}})));

    replaceable Battery.Cells.Aging.Variants.NoAging agingModel constrainedby
      Battery.Cells.Aging.Partials.PartialAgingModel(
      final SOH_init=SOH_init, final SOHR_init=SOHR_init) "Aging model in cell"
      annotation (
      Dialog(group="Replaceable Models"),
      choicesAllMatching=true,
      Placement(transformation(extent={{-20,-80},{20,-40}})));

    parameter Boolean assertMinSOC=false "Terminate simulation if SOC is < SOC_minAssert" annotation (
      Evaluate=true,
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(group="SOC Options", tab="Assert options"));
    parameter Real SOC_minAssert=-0.1 "Minimal value for SOC" annotation (Dialog(
        enable=assertMinSOC,
        group="SOC Options",
        tab="Assert options"));

    parameter Boolean assertMaxSOC=false "Terminate simulation if SOC is > SOC_maxAssert " annotation (
      Evaluate=true,
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(group="SOC Options", tab="Assert options"));
    parameter Real SOC_maxAssert=1.1 "Maximal value for SOC" annotation (Dialog(
        enable=assertMaxSOC,
        group="SOC Options",
        tab="Assert options"));

    parameter Boolean assertMinV=false "Terminate simulation if v is < V_minAssert" annotation (
      Evaluate=true,
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(group="Voltage Options", tab="Assert options"));
    parameter Real V_minAssert=2.7 "Minimal value for v" annotation (Dialog(
        enable=assertMinV,
        group="Voltage Options",
        tab="Assert options"));

    parameter Boolean assertMaxV=false "Terminate simulation if v is > V_maxAssert " annotation (
      Evaluate=true,
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(group="Voltage Options", tab="Assert options"));
    parameter Real V_maxAssert=4.3 "Maximal value for v" annotation (Dialog(
        enable=assertMaxV,
        group="Voltage Options",
        tab="Assert options"));

    parameter Boolean assertMinI=false "Terminate simulation if i is < I_minAssert" annotation (
      Evaluate=true,
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(group="Current Options", tab="Assert options"));
    parameter Real I_minAssert=-15*electricModel.C_nominal/3600 "Minimal value for Current (current sign is negative when discharging)"
      annotation (Dialog(
        enable=assertMinI,
        group="Current Options",
        tab="Assert options"));

    parameter Boolean assertMaxI=false "Terminate simulation if i is > I_maxAssert (current sign is positive when charging)"
      annotation (
      Evaluate=true,
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(group="Current Options", tab="Assert options"));
    parameter Real I_maxAssert=4*electricModel.C_nominal/3600 "Maximal value for Current" annotation (Dialog(
        enable=assertMaxI,
        group="Current Options",
        tab="Assert options"));

    parameter Boolean assertMinT=false "Terminate simulation if T is < T_minAssert" annotation (
      Evaluate=true,
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(group="Temperature Options", tab="Assert options"));
    parameter Modelica.Units.SI.Temperature T_minAssert=273.14 - 35 "Minimal value for T"
      annotation (Dialog(
        enable=assertMinT,
        group="Temperature Options",
        tab="Assert options"));

    parameter Boolean assertMaxT=false "Terminate simulation if T is > T_maxAssert " annotation (
      Evaluate=true,
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(group="Temperature Options", tab="Assert options"));
    parameter Modelica.Units.SI.Temperature T_maxAssert=273.14 + 65 "Maximal value for T"
      annotation (Dialog(
        enable=assertMaxT,
        group="Temperature Options",
        tab="Assert options"));

    parameter Boolean assertMinSOH=false "Terminate simulation if SOH is < SOH_minAssert" annotation (
      Evaluate=true,
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(group="SOH Options", tab="Assert options"));
    parameter Real SOH_minAssert=0.7 "Minimal value for SOH" annotation (Dialog(
        enable=assertMinSOH,
        group="SOH Options",
        tab="Assert options"));

    parameter Boolean assertMaxSOH=false "Terminate simulation if SOH is > SOH_maxAssert " annotation (
      Evaluate=true,
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(group="SOH Options", tab="Assert options"));
    parameter Real SOH_maxAssert=1 "Maximal value for SOH" annotation (Dialog(
        enable=assertMaxSOH,
        group="SOH Options",
        tab="Assert options"));

  equation

    if assertMinSOC then
      assert(SOC_minAssert <= electricModel.SOC, "SOC below lower boundary of " + String(SOC_minAssert));
    end if;
    if assertMinT then
      assert(T_minAssert <= electricModel.T, "Temperature below lower boundary of " + String(T_minAssert) + " K");
    end if;
    if assertMinI then
      assert(I_minAssert <= electricModel.i, "Current below lower boundary of " + String(I_minAssert) + " A");
    end if;
    if assertMinV then
      assert(V_minAssert <= electricModel.v, "Voltage below lower boundary of " + String(V_minAssert) + " V");
    end if;
    if assertMinSOH then
      assert(SOH_minAssert <= electricModel.SOH, "SOH below lower boundary of " + String(SOH_minAssert));
    end if;
    if assertMaxSOC then
      assert(SOC_maxAssert >= electricModel.SOC, "SOC above upper boundary of " + String(SOC_maxAssert));
    end if;
    if assertMaxT then
      assert(T_maxAssert >= electricModel.T, "Temperature above upper boundary of " + String(T_maxAssert) + " K");
    end if;
    if assertMaxI then
      assert(I_maxAssert >= electricModel.i, "Current above upper boundary of " + String(I_maxAssert) + " A");
    end if;
    if assertMaxV then
      assert(V_maxAssert >= electricModel.v, "Voltage above upper boundary of " + String(V_maxAssert) + " V");
    end if;
    if assertMaxSOH then
      assert(SOH_maxAssert >= electricModel.SOH, "SOH above upper boundary of " + String(SOH_maxAssert));
    end if;

    connect(p, electricModel.p) annotation (Line(
        points={{-100,0},{-20,0}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(electricModel.n, n) annotation (Line(
        points={{20,0},{100,0}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(cellBus, agingModel.cellBus) annotation (Line(
        points={{0,-100},{0,-80}},
        color={255,204,51},
        thickness=0.5));
    connect(electricModel.cellBus, agingModel.cellBus) annotation (Line(
        points={{0,-20},{0,-30},{40,-30},{40,-90},{0,-90},{0,-80}},
        color={255,204,51},
        thickness=0.5));
    connect(thermalModel.cellBus, agingModel.cellBus) annotation (Line(
        points={{0,40},{0,30},{40,30},{40,-90},{0,-90},{0,-80}},
        color={255,204,51},
        thickness=0.5));
    connect(thermalModel.cylindricHeatPort, heatPort) annotation (Line(points={{0,80},{0,100}}, color={191,0,0}));
    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}})),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}}), graphics={
          Rectangle(
            extent={{-10,72},{10,42}},
            lineColor={255,255,255},
            fillColor={0,86,134},
            fillPattern=FillPattern.Solid,
            radius=1),
          Rectangle(
            extent={{-40,54},{40,-62}},
            lineColor={255,255,255},
            fillColor={0,86,134},
            fillPattern=FillPattern.Solid,
            radius=2)}),
      Documentation(info="<html>
<p>The CellTemplate class is the template to extend from in order to create new battery cell models. For more information on the structure of the cell model and its submodels see the <a href=\"modelica://Battery.UsersGuide.Cells\">cell</a> section of the user guide. </p>
<p>The cell model consists of three replaceable sub models described in the list below: </p>
<ol>
<li>A thermal model representing the thermal behavior of the battery cell. The package and all the information about the thermal models can be found <a href=\"modelica://Battery.Cells.Thermal\">here</a>. </li>
<li>An electric model representing the electrical behavior of the battery cell. The package and the information about the electrical models can be found <a href=\"modelica://Battery.Cells.Electric\">here</a>. </li>
<li>An aging model representing the aging behavior of the battery cell. The package and all the information about the aging models can be found <a href=\"modelica://Battery.Cells.Aging\">here</a>. </li>
</ol>
<p>The parameters of a cell model based on this template are the simulation start values for temperature <i>T</i> and state of charge <i>SOC</i>, the number of discretized elements of the thermal model <i>N_verticalElements</i> in z-direction (height), the number of surface heat ports on the perimeter of the cell <i>N_surface</i> and the assert options defining the range of operation (See the tab <i>Assert options</i> in the parmeter dialog). The parameter <span style=\"font-family: Courier New;\">N_verticalElements</span> allows to discretize the volume of the cell core in the thermal model. </p>
<p><b>Geometry and orientation</b> </p>
<p>All geometrical parameters are defined in the thermal model of the cell and only influence the thermal model with its themal capacity and thermal resistances, none of these parameters influence the electrical model neither for the ohmic resistance nor for the capacity. </p>
<p>The geomtry of a cylindric cell is shown in <b>Fig. 1</b>. </p>
<p><b>Fig. 1: </b>Coordinate system and dimensions for cylindric cell models </p>
<table cellspacing=\"0\" cellpadding=\"2\" border=\"0\"><tr>
<td><p><img src=\"modelica://Battery/Resources/Images/UsersGuide/Cells/CellCylindric.png\" alt=\"coordinate_system_cell.png\"/> </p></td>
</tr>
</table>
<p><br><h4>Thermal connection</h4></p>
<p>Each battery cell features a single connector for all available heat ports. The model of this connector is called <a href=\"Battery.Common.Interfaces.CylindricHeatPort\">CylindricHeatPort</a>. </p>
</html>"));
  end CylindricBat;
  annotation (uses(Modelica(version="4.0.0"), Battery(version="2.2.0")));
end ETBATTERY;
