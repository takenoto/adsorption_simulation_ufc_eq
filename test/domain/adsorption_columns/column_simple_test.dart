import 'package:adsorption_column_sim_renan_flutter/domain/domain.dart';
import 'package:flutter_test/flutter_test.dart';

void main(){

  test('can run column', (){
    final s1Key = 'substance_1';
    final s2Key = 'substance_2';
    final L = 1; //cm
    final divisions = 5;
    final spaceVector = UnidimensionalDiscretization(List.generate(divisions, (index) => L*(index/divisions)));
    final zeros = UnidimensionalDiscretization(spaceVector.clone().values.map((e) => 0.0).toList());
    final ones = UnidimensionalDiscretization(spaceVector.clone().values.map((e) => 1.0).toList());
    final concsZero = UnidimensionalConcentration<String>(
      {
        s1Key: ones,
        s2Key: zeros,
      }
    );
    final initialState = UnidimensionalColumnState(fluidPhaseConcs: concsZero, adsorbedPhaseConcs: concsZero, space: spaceVector);
    final idealColumn = IdealUnidimensionalChromatographicColumn(
      initialState: initialState,
      u: 0.5,
      nu: 10,
      porosity: 0.3,
      dax: 0.5,
      dqDt: LumpedDqDtCalculator<String>(
        {
          s1Key: ({required double C, required double q}) => C*0.1/100,
          s2Key: ({required double C, required double q}) => C*0.7/220,
        }
      )
    );

    final derivatives = idealColumn.calculateDerivatives(initialState);
    // print(derivatives.derivativesC[s1Key]!);
    // print(derivatives.derivativesC[s2Key]!);
    // print(derivatives.derivativesQ[s1Key]!);
    // print(derivatives.derivativesQ[s2Key]!);
    
    //FIXME problema óbvio detectado: a concentração está aumentando tanto na fase líquida quanto na adsorvida ?????
    final newState = idealColumn.newStateUsingDerivatives(current: initialState, derivatives: derivatives, dtInSec: 100);
    print(initialState.fluidPhaseConcs.concs[s1Key]!.values);
    print(newState.fluidPhaseConcs.concs[s1Key]!.values);
    print(initialState.adsorbedPhaseConcs.concs[s1Key]!.values);
    print(newState.adsorbedPhaseConcs.concs[s1Key]!.values); 
    
  });
}