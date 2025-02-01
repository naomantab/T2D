// Grab unique values for each column
function getUniqueValuesFromColumn() {
    let unique_col_values_dict = {};

    let allFilters = document.querySelectorAll(".table-filter");
    allFilters.forEach((filter_i) => {
        let col_index = filter_i.parentElement.getAttribute("col-index");

        const rows = document.querySelectorAll("#snp-table > tbody > tr");

        rows.forEach((row) => {
            let cell_value = row.querySelector("td:nth-child(" + col_index + ")").innerText.trim();

            if (col_index in unique_col_values_dict) {
                if (!unique_col_values_dict[col_index].includes(cell_value)) {
                    unique_col_values_dict[col_index].push(cell_value);
                }
            } else {
                unique_col_values_dict[col_index] = [cell_value];
            }
        });
    });

    updateSelectOptions(unique_col_values_dict);
}

// Add <option> tag to each filter dropdown based on unique values
function updateSelectOptions(unique_col_values_dict) {
    let allFilters = document.querySelectorAll(".table-filter");

    allFilters.forEach((filter_i) => {
        let col_index = filter_i.parentElement.getAttribute("col-index");
        filter_i.innerHTML = '<option value="all">All</option>'; // Reset filter options

        unique_col_values_dict[col_index].forEach((value) => {
            filter_i.innerHTML += `<option value="${value}">${value}</option>`;
        });
    });
}

// Filter table rows based on selected filter values
function filter_rows() {
    let allFilters = document.querySelectorAll(".table-filter");
    let filter_value_dict = {};

    allFilters.forEach((filter_i) => {
        let col_index = filter_i.parentElement.getAttribute("col-index");
        let value = filter_i.value;

        if (value !== "all") {
            filter_value_dict[col_index] = value;
        }
    });

    const rows = document.querySelectorAll("#snp-table tbody tr");
    
    rows.forEach((row) => {
        let display_row = true;

        for (let col_i in filter_value_dict) {
            let row_cell_value = row.querySelector("td:nth-child(" + col_i + ")").innerText.trim();
            let filter_value = filter_value_dict[col_i];

            if (row_cell_value.toLowerCase().indexOf(filter_value.toLowerCase()) === -1) {
                display_row = false;
                break;
            }
        }

        row.style.display = display_row ? "table-row" : "none";
    });
}

// Event Listener for filter changes
document.addEventListener("DOMContentLoaded", () => {
    getUniqueValuesFromColumn();
    document.querySelectorAll(".table-filter").forEach((filter) => {
        filter.addEventListener("change", filter_rows);
    });
});
