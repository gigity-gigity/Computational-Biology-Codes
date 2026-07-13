# import numpy as np
# import matplotlib.pyplot as plt
#
# carbon=['C32','C33','C34','C35','C36','C37','C38','C39','C310','C311','C312','C313','C314','C315','C316']
# order_parameter={}
# with open("/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/300_sn1.orderparSCD") as handle:
#     for carbon_no in carbon:
#         average = []
#         handle.seek(0)
#         for lines in handle:
#             line=lines.strip()
#             if carbon_no in line:
#                 average.append(float(line.split('\t')[-1]))
#         if len(average) == 0:
#             print(f"No data found for {carbon_no}")
#             continue
#         mean_value = round(sum(average) / len(average),3)
#         variance = sum([((x - mean_value) ** 2) for x in average]) / len(average)
#         std_deviation = round(variance ** 0.5,3)
#         order_parameter[carbon_no]=mean_value,std_deviation
# print(order_parameter)
#
# carbon_numbers, values = zip(*order_parameter.items())
# means, std_devs = zip(*values)
# # Plotting
# fig, ax = plt.subplots()
# fig.patch.set_facecolor('white')  # Set the background color of the entire plot to white
# ax.plot(carbon_numbers, means, marker='o', color='blue', label='Mean Values')
# ax.errorbar(carbon_numbers, means, yerr=std_devs, fmt='none', ecolor='gray', capsize=5, label='Standard Deviation')
#
# # Adding labels and title
# ax.set_xlabel('Carbon Numbers')
# ax.set_ylabel('SCD')
# ax.set_title('Order parameter of sn-1')
# # Set y-axis ticks with a 0.05 interval
# ax.set_yticks(np.arange(-0.2, 0.5, 0.05))
# ax.set_yticklabels([round(i, 2) for i in ax.get_yticks()])  # Adjust labels for better readability
#
# ax.set_xticks(carbon_numbers)
# ax.set_xticklabels(carbon_numbers, rotation=90, ha='right')  # Adjust rotation for better readability
# ax.legend()
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
#
# plt.savefig("sn1.jpeg", bbox_inches="tight")
# plt.show()
# plt.plot()
#####################################################################################################################
                                 ############### SN-2 plot #############
import numpy as np
import matplotlib.pyplot as plt
carbon=['C22','C23','C24','C25','C26','C27','C28','C29','C210','C211','C212','C213','C214','C215','C216','C217','C218']
order_parameter={}
with open("/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/300_sn2.orderparSCD") as handle:
    for carbon_no in carbon:
        average = []
        handle.seek(0)
        for lines in handle:
            line=lines.strip()
            if carbon_no in line:
                average.append(float(line.split('\t')[-1]))
        if len(average) == 0:
            print(f"No data found for {carbon_no}")
            continue
        mean_value = round(sum(average) / len(average),3)
        variance = sum([((x - mean_value) ** 2) for x in average]) / len(average)
        std_deviation = round(variance ** 0.5,3)
        order_parameter[carbon_no]=mean_value,std_deviation
print(order_parameter)

carbon_numbers, values = zip(*order_parameter.items())
means, std_devs = zip(*values)
# Plotting
fig, ax = plt.subplots()
fig.patch.set_facecolor('white')  # Set the background color of the entire plot to white
ax.plot(carbon_numbers, means, marker='o', color='blue', label='Mean Values')
ax.errorbar(carbon_numbers, means, yerr=std_devs, fmt='none', ecolor='gray', capsize=5, label='Standard Deviation')

# Adding labels and title
ax.set_xlabel('Carbon Numbers')
ax.set_ylabel('SCD')
ax.set_title('Order parameter of sn-2')
# Set y-axis ticks with a 0.05 interval
ax.set_yticks(np.arange(-0.4, 0.6, 0.05))
ax.set_yticklabels([round(i, 2) for i in ax.get_yticks()])  # Adjust labels for better readability

ax.set_xticks(carbon_numbers)
ax.set_xticklabels(carbon_numbers, rotation=90, ha='right')  # Adjust rotation for better readability
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig("sn2.jpeg", bbox_inches="tight")
plt.show()
plt.plot()






