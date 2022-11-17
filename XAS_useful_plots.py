#Delay vs I1_off, I0_1_off and I0_1_off, rebinned
plt.plot( Delays_off_time_rebinned, I1_off_time_rebinned)
plt.plot( Delays_off_time_rebinned, I0_1_off_time_rebinned)
plt.plot( Delays_off_time_rebinned, I0_2_off_time_rebinned)
plt.legend(['I1', 'I0_1', 'I0_2'])

#Delay vs standard deviations of I1_off, I0_1_off and I0_1_off, rebinned
plt.plot( Delays_off_time_rebinned, I1_off_time_rebinned_std)
plt.plot( Delays_off_time_rebinned, I0_1_off_time_rebinned_std)
plt.plot( Delays_off_time_rebinned, I0_2_off_time_rebinned_std)
plt.legend(['I1', 'I0_1', 'I0_2'])

#Delay vs I1_off, I0_1_off and I0_1_off with error bars, rebinned
plt.errorbar( Delays_off_time_rebinned, I1_off_time_rebinned, yerr=I1_off_time_rebinned_std, uplims=True, lolims=True)
plt.errorbar( Delays_off_time_rebinned, I0_1_off_time_rebinned, yerr=I0_1_off_time_rebinned_std, uplims=True, lolims=True)
plt.errorbar( Delays_off_time_rebinned, I0_2_off_time_rebinned, yerr=I0_2_off_time_rebinned_std, uplims=True, lolims=True)
plt.legend(['I1', 'I0_1', 'I0_2'])

#Delay vs I1_on, I0_1_on and I0_1_on, rebinned
plt.plot( Delays_on_time_rebinned, I1_on_time_rebinned)
plt.plot( Delays_on_time_rebinned, I0_1_on_time_rebinned)
plt.plot( Delays_on_time_rebinned, I0_2_on_time_rebinned)
plt.legend(['I1', 'I0_1', 'I0_2'])

#Delay vs standard deviations of I1_on, I0_1_on and I0_1_on, rebinned
plt.plot( Delays_on_time_rebinned, I1_on_time_rebinned_std)
plt.plot( Delays_on_time_rebinned, I0_1_on_time_rebinned_std)
plt.plot( Delays_on_time_rebinned, I0_2_on_time_rebinned_std)
plt.legend(['I1', 'I0_1', 'I0_2'])

#Delay vs I1_on, I0_1_on and I0_1_on with error bars, rebinned
plt.errorbar( Delays_on_time_rebinned, I1_on_time_rebinned, yerr=I1_on_time_rebinned_std, uplims=True, lolims=True)
plt.errorbar( Delays_on_time_rebinned, I0_1_on_time_rebinned, yerr=I0_1_on_time_rebinned_std, uplims=True, lolims=True)
plt.errorbar( Delays_on_time_rebinned, I0_2_on_time_rebinned, yerr=I0_2_on_time_rebinned_std, uplims=True, lolims=True)
plt.legend(['I1', 'I0_1', 'I0_2'])

#Delay vs normalized I1_on, I1_off and their difference, rebinned
plt.plot( Delays_on_time_rebinned, (I1_on_time_rebinned/(I0_1_on_time_rebinned+I0_2_on_time_rebinned)) )
plt.plot( Delays_off_time_rebinned, (I1_off_time_rebinned/(I0_1_off_time_rebinned+I0_2_off_time_rebinned)) )
plt.plot( Delays_on_time_rebinned[0:85], (I1_on_time_rebinned[0:85]/(I0_1_on_time_rebinned[0:85]+I0_2_on_time_rebinned[0:85]))-(I1_off_time_rebinned[0:85]/(I0_1_off_time_rebinned[0:85]+I0_2_off_time_rebinned[0:85])) )

#Delay vs normalized I1_on with error bars, error propagation rules applied, rebinned
I1 = I1_on_time_rebinned/(I0_1_on_time_rebinned+I0_2_on_time_rebinned)
I1_std = np.sqrt( (I1_on_time_rebinned_std/(I0_1_on_time_rebinned+I0_2_on_time_rebinned))**2 + (-I1_on_time_rebinned*I0_1_on_time_rebinned_std/((I0_1_on_time_rebinned+I0_2_on_time_rebinned)**2))**2 + (-I1_on_time_rebinned*I0_2_on_time_rebinned_std/((I0_1_on_time_rebinned+I0_2_on_time_rebinned)**2))**2 )
plt.errorbar( Delays_on_time_rebinned, I1, yerr=I1_std, uplims=True, lolims=True)

#I1_on vs I0_1_on, I0_2_on
plt.plot( I1_on_time_rebinned, I0_1_on_time_rebinned,'.')
plt.plot( I1_on_time_rebinned, I0_2_on_time_rebinned,'.')
