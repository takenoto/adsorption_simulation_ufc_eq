import 'package:adsorption_column_sim_renan_flutter/domain/discretization/discretization.dart';

import 'calculators.dart';

class UnidimensionalConcentration<T> {
  ///A key é a substância de interesse
  ///Pode ser uma String com o casNumber, um int que a represente, etc.
  final Map<T, UnidimensionalDiscretization> concs;
  const UnidimensionalConcentration(this.concs);

  UnidimensionalConcentration clone() {
    final newMap = <T, UnidimensionalDiscretization>{};
    concs.forEach((key, value) {
      newMap[key] = value.clone();
    });
    return UnidimensionalConcentration<T>(newMap);
  }
}

///The length of the lists of each fluid and adsorbed phase must be equals to the [space] length.
class UnidimensionalColumnState {
  final UnidimensionalConcentration<String> fluidPhaseConcs;
  final UnidimensionalConcentration<String> adsorbedPhaseConcs;
  final UnidimensionalDiscretization space;
  UnidimensionalColumnState(
      {required this.fluidPhaseConcs,
      required this.adsorbedPhaseConcs,
      required this.space});

  UnidimensionalColumnState clone(){
      return  UnidimensionalColumnState(
        space: space.clone(),
        adsorbedPhaseConcs: adsorbedPhaseConcs.clone()
            as UnidimensionalConcentration<String>,
        fluidPhaseConcs: fluidPhaseConcs.clone()
            as UnidimensionalConcentration<String>);
  }
}

class IdealUnidimChromColumnDerivatives{
  final Map<String, List<double>> derivativesC;
  final Map<String, List<double>> derivativesQ;
  IdealUnidimChromColumnDerivatives(this.derivativesC, this.derivativesQ);
}

//PENSAR: Como seria a coluna ideal simples?
class IdealUnidimensionalChromatographicColumn {
  //Os valores iniciais para cada substância/coisa
  ///Valores da fase fluido
  final UnidimensionalColumnState initialState;
  final double dax;
  final double porosity;
  final double u;
  final double nu;

  final _dcDt = UnidimensionalDcDt();
  final singleDerivativeCalc = SingleDerivativeCalculatorFiniteDifference();
  final secondDerivativeCalc = SecondDerivativeCalculatorFiniteDifference();
  final LumpedDqDtCalculator dqDt;

  IdealUnidimensionalChromatographicColumn(
      {required this.initialState,
      required this.dqDt,
      required this.dax,
      required this.porosity,
      required this.u,
      required this.nu});

  ///Os valores opcionais podem ser usados para  sobrescrever os valores passados no criador dessa coluna.
  ///[dt] é em segundos
  IdealUnidimChromColumnDerivatives calculateDerivatives(UnidimensionalColumnState from,
      {double? porosity, double? u, double? dax, double? nu}) {
        //Prepara um novo estado, limpo, para usar.
    UnidimensionalColumnState newState = from.clone();

    //Determina quais as keys existentes para as substâncias
    final substancesKeys = newState.fluidPhaseConcs.concs.keys.toList();
    //As derivadas de concentração, para cada substância (separadas por chave) no leito e no adsorvente
    final derivativesC = <String, List<double>>{};
    final derivativesQ = <String, List<double>>{};
    //Inicializa com vetores vazios
    for (var s in substancesKeys) {
      derivativesC[s] = <double>[];
      derivativesQ[s] = <double>[];
    }

    final space = newState.space;
    //Passa por todos os valores do step anterior, e descobre como estariam no step seguinte
    for (int i = 0; i < space.values.length; i++) {
      //Espaço local:
      double? previousZ;
      if (i > 0) previousZ = space.values[i - 1];
      double z = space.values[i];
      double? nextZ;
      if (i < space.values.length-1) nextZ = space.values[i + 1];
      //Faz os cálculos para cada substância, uma a uma
      for (var substanceKey in substancesKeys) {
        //Para a posição atual:
        //1) Checa se é uma condição de contorno e, portanto, não muda
        //TODO break;
        final concentrationProfile =
            from.fluidPhaseConcs.concs[substanceKey]!.values;
        //2)Se não for condição de contorno (ou seja, se variar) faz o cálculo do novo valor
        //Calcula o valor de d²c/dx²
        //Se for em qualquer ponto interno, faz derivada central
        //Se for em qualquer outro, faz a que for possível
        double? d2c_dx2, dc_dx, dq_dt;
        final dq_dtFunc = dqDt.forSubstance(substanceKey);
        //TODO checa se nessa posição d²c/dx² é condição de contorno. for, usa o valor. Se não, faz os cálculos.
        if (previousZ != null && nextZ != null) {
          //d²c/dx²
          //É possível usar os dois para calcular por diferenças centrais
          final s = SecondDerivativeInput(
              f_minus_1: concentrationProfile[i - 1],
              f_0: concentrationProfile[i],
              f_plus_1: concentrationProfile[i + 1],
              h: (nextZ - previousZ) / 2);
          d2c_dx2 = secondDerivativeCalc.d2ydx2(s);
        } 
        
        if (previousZ != null) {
          //O anterior é zero, então precisa usar o próximo
          //dc/dx
          dc_dx = singleDerivativeCalc.difference(
              ValuePair(concentrationProfile[i - 1], concentrationProfile[i]),
              ValuePair(previousZ, z));
        } else if (nextZ != null) {
          //dc/dx
          dc_dx = singleDerivativeCalc.difference(
              ValuePair(concentrationProfile[i], concentrationProfile[i + 1]),
              ValuePair(z, nextZ));
        }

        dq_dt = dq_dtFunc(
            C: concentrationProfile[i],
            q: from.adsorbedPhaseConcs.concs[substanceKey]!.values[i]);

        //2.1) Calcula o valor de dcdt:
        final dcDt = _dcDt.dc_dt(
            dax: dax ?? this.dax,
            d2c_dx2: d2c_dx2??0,
            nu: nu ?? this.nu,
            u: u ?? this.u,
            dc_dx: dc_dx!,
            porosity: porosity ?? this.porosity,
            dq_dt: dq_dt);

        //Adiciona o valor de dcDt das substâncias da fase móvels
        derivativesC[substanceKey]?.add(dcDt);
        //Adiciona o valor de dQ/dt das substâncias na fase
        derivativesQ[substanceKey]?.add(dq_dt);
      }
    }
    return IdealUnidimChromColumnDerivatives(
      derivativesC,
      derivativesQ
    );
    
  }
  
  //[dt] must be IN SECONDS!
  UnidimensionalColumnState newStateUsingDerivatives(
    {required UnidimensionalColumnState current,
    required IdealUnidimChromColumnDerivatives derivatives,
    required double dtInSec}
  ){
    final dt = dtInSec;
    final newState = current.clone();
    
    for(int i=0; i<newState.space.values.length; i++){
      for(var key in newState.fluidPhaseConcs.concs.keys){
         newState.fluidPhaseConcs.concs[key]!.values[i] += derivatives.derivativesC[key]![i]*dt;
         newState.adsorbedPhaseConcs.concs[key]!.values[i] += derivatives.derivativesQ[key]![i]*dt;
      }
    }
    return newState;
  }
}

//------------------------------------------------
//Daqui pra baixo ignore acho que não prestou

// abstract class Adsorbent {
//   final String name;
//   final String? casNumber;
//   Adsorbent(this.name, [this.casNumber]);

//   double getQAtEquilibriumForSubstance(String substanceCasNumber,
//       double temperatureInKelvin, double pressureInPa) {
//     return 0;
//   }
// }

// //Acho que nem vou usar isso...
// class AdsorbentParticle {
//   final Adsorbent adsorbent;
//   AdsorbentParticle(this.adsorbent);
// }

// class ColumnRegion {
//   double temperatureInKelvin;
//   double pressureInPa;
//   ColumnRegion(this.temperatureInKelvin, this.pressureInPa);
// }

// enum ColumnSubstanceLocation { fluid, adsorbent }

// class CircularColumnRegion extends ColumnRegion {
//   Map<String, double> fluidSubstancesConc;
//   Map<String, double> adsorbentSubstancesConc;
//   final Adsorbent adsorbent;
//   CircularColumnRegion(
//       {required this.fluidSubstancesConc,
//       required this.adsorbentSubstancesConc,
//       required this.adsorbent,
//       required double temperatureInKelvin,
//       required double pressureInPa})
//       : super(temperatureInKelvin, pressureInPa);

//   double getSubstanceConcentrationAt(
//       String substanceCas, ColumnSubstanceLocation at) {
//     double? substanceConc = 0;
//     switch (at) {
//       case ColumnSubstanceLocation.fluid:
//         substanceConc = fluidSubstancesConc[substanceCas];
//         break;
//       case ColumnSubstanceLocation.adsorbent:
//         substanceConc = adsorbentSubstancesConc[substanceCas];
//         break;
//       default:
//         break;
//     }
//     //Por padrão, se não encontrar nada, é zero.
//     return substanceConc ?? 0;
//   }

//   ///Temperature must be in Kelvin
//   ///Pressure must be in Pa
//   double getQAtEquilibriumForSubstance(
//       String substanceCasNumber, double t, double p) {
//     return adsorbent.getQAtEquilibriumForSubstance(substanceCasNumber, t, p);
//   }
// }

// class CircularColumn {
//   CircularColumn._(
//       {required this.dz, required this.radiusInCm, required this.regions});

//   final List<CircularColumnRegion> regions;
//   final double radiusInCm;
//   final double dz;

//   //Se quiser uma nova coluna (snapshot), implementar um copyWith
//   //Retorna o valor de dc/dt para cada substância nesse delta de espaço.
//   //A chave do mapa é a id da substância.
//   //O valor é uma lista com seu dc/dt em cada ponto do espaço.

//   //Eu acho que sses valores deviam ser passados é localmente zzz

//   //Talvez eu já pudesse editar em algo nesse estilo, e criar helpers
//   //Pra tornar mais fácil apenas.
//   //A resposta é nesse estilo:
//   ///{
//   // space:{
//   // 'z': [] //lista,
//   // 'x': null
//   // 'y': null
//   // 'radius': null
//   // 'angular': null
//   //},
//   //concentration:{
//   // [
//   // {
//   // "substance_name":  'água',
//   //  'cas_number': 'NUMERO CAS 123123213',
//   // 'z': [] //lista com os valores e assim por diante
//   // }
//   //
//   //]
//   //}
//   //
//   //}
//   Map<String, dynamic> step(double dtInMs,
//       {required double dax,
//       required double d2c_dx2,
//       required double u,
//       required double dc_dx,
//       required double porosidade,
//       required double dq_dt}) {
//     final profiles = <String, dynamic>{};
//     //1º Para cada região:
//     for (int i = 0; i < regions.length; i++) {
//       //Calcula a adsorção/dessorção de todas as substâncias de interesse
//       //TODO calcula o dcDt mesmo e deixa pra integrar depois.
//       //Seguindo o modelo mais simples dos slides do prof. Amaro. Ainda não entendi bem qual é esse...
//       //Calcula o dcDt
//       var dcDt =
//           dax * d2c_dx2 - u * dc_dx + ((1 - porosidade) / porosidade) * dq_dt;

//       //2º e interage com a próxima região através de
//     }

//     return profiles;
//   }

//   factory CircularColumn.create(
//       {required int numberOfRegions,
//       required double radiusInCm,
//       required double columnLength,
//       required Adsorbent adsorbent,
//       double initialTemperatureInKelvin = 293,
//       double initialPressureInPa = 1e5}) {
//     //Descobre o número de regiões necessárias para representar essa coluna no dz desejado
//     //Gera as regiões, em ordem, do topo para o fundo:
//     final regions = List.generate(numberOfRegions, (i) {
//       return CircularColumnRegion(
//           fluidSubstancesConc: {},
//           adsorbentSubstancesConc: {},
//           adsorbent: adsorbent,
//           temperatureInKelvin: initialTemperatureInKelvin,
//           pressureInPa: initialPressureInPa);
//     });
//     return CircularColumn._(
//         dz: columnLength / numberOfRegions,
//         radiusInCm: radiusInCm,
//         regions: regions);
//   }
// }
